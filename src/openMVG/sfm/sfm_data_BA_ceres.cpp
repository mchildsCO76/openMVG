// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;

/// Create the appropriate cost functor according the provided input camera intrinsic model.
/// The residual can be weighetd if desired (default 0.0 means no weight).
ceres::CostFunction * IntrinsicsToCostFunction
(
  IntrinsicBase * intrinsic,
  const Vec2 & observation,
  const double weight
)
{
  switch(intrinsic->getType())
  {
    case PINHOLE_CAMERA:
        return ResidualErrorFunctor_Pinhole_Intrinsic::Create(observation, weight);
     break;
    case PINHOLE_CAMERA_RADIAL1:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K1::Create(observation, weight);
    break;
    case PINHOLE_CAMERA_RADIAL3:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Radial_K3::Create(observation, weight);
    break;
    case PINHOLE_CAMERA_BROWN:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Brown_T2::Create(observation, weight);
    break;
    case PINHOLE_CAMERA_FISHEYE:
      return ResidualErrorFunctor_Pinhole_Intrinsic_Fisheye::Create(observation, weight);
    default:
      return nullptr;
  }
}

Bundle_Adjustment_Ceres::BA_Ceres_options::BA_Ceres_options
(
  const bool bVerbose,
  bool bmultithreaded
)
: bVerbose_(bVerbose),
  nb_threads_(1),
  parameter_tolerance_(1e-8), //~= numeric_limits<float>::epsilon()
  bUse_loss_function_(true)
{
  #ifdef OPENMVG_USE_OPENMP
    nb_threads_ = omp_get_max_threads();
  #endif // OPENMVG_USE_OPENMP
  if (!bmultithreaded)
    nb_threads_ = 1;

  bCeres_summary_ = false;

  // Default configuration use a DENSE representation
  linear_solver_type_ = ceres::DENSE_SCHUR;
  preconditioner_type_ = ceres::JACOBI;
  // If Sparse linear solver are available
  // Descending priority order by efficiency (SUITE_SPARSE > CX_SPARSE > EIGEN_SPARSE)
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
  {
    sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
    linear_solver_type_ = ceres::SPARSE_SCHUR;
  }
  else
  {
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
    {
      sparse_linear_algebra_library_type_ = ceres::CX_SPARSE;
      linear_solver_type_ = ceres::SPARSE_SCHUR;
    }
    else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
    {
      sparse_linear_algebra_library_type_ = ceres::EIGEN_SPARSE;
      linear_solver_type_ = ceres::SPARSE_SCHUR;
    }
  }
}


Bundle_Adjustment_Ceres::Bundle_Adjustment_Ceres
(
  Bundle_Adjustment_Ceres::BA_Ceres_options options
)
: ceres_options_(options)
{}

Bundle_Adjustment_Ceres::BA_Ceres_options &
Bundle_Adjustment_Ceres::ceres_options()
{
  return ceres_options_;
}

bool Bundle_Adjustment_Ceres::Adjust
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
 for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
  {
    const IndexT indexPose = itPose->first;

    const Pose3 & pose = itPose->second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
    itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
  {
    const IndexT indexCam = itIntrinsic->first;

    if (isValid(itIntrinsic->second->getType()))
    {
      map_intrinsics[indexCam] = itIntrinsic->second->getParams();

      double * parameter_block = &map_intrinsics[indexCam][0];
      problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
      if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
      {
        // set the whole parameter block as constant for best performance
        problem.SetParameterBlockConstant(parameter_block);
      }
      else
      {
        const std::vector<int> vec_constant_intrinsic =
          itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
        if (!vec_constant_intrinsic.empty())
        {
          ceres::SubsetParameterization *subset_parameterization =
            new ceres::SubsetParameterization(
              map_intrinsics[indexCam].size(), vec_constant_intrinsic);
          problem.SetParameterization(parameter_block, subset_parameterization);
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks!= sfm_data.structure.end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;

    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(itObs->first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

      if (cost_function)
        problem.AddResidualBlock(cost_function,
          p_LossFunction,
          &map_intrinsics[view->id_intrinsic][0],
          &map_poses[view->id_pose][0],
          iterTracks->second.X.data());
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(iterTracks->second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (Landmarks::iterator iterGCPTracks = sfm_data.control_points.begin();
      iterGCPTracks!= sfm_data.control_points.end(); ++iterGCPTracks)
    {
      const Observations & obs = iterGCPTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            itObs->second.x,
            options.control_point_opt.weight);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            nullptr,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            iterGCPTracks->second.X.data());
      }
      // Set the 3D point as FIXED (it's a GCP)
      problem.SetParameterBlockConstant(iterGCPTracks->second.X.data());
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
  ceres::Solver::Options ceres_config_options;
  ceres_config_options.preconditioner_type = ceres_options_.preconditioner_type_;
  ceres_config_options.linear_solver_type = ceres_options_.linear_solver_type_;
  ceres_config_options.sparse_linear_algebra_library_type = ceres_options_.sparse_linear_algebra_library_type_;
  ceres_config_options.minimizer_progress_to_stdout = false;
  ceres_config_options.logging_type = ceres::SILENT;
  ceres_config_options.num_threads = ceres_options_.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options_.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options_.parameter_tolerance_;

  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options_.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options_.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (Poses::iterator itPose = sfm_data.poses.begin();
        itPose != sfm_data.poses.end(); ++itPose)
      {
        const IndexT indexPose = itPose->first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = itPose->second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (Intrinsics::iterator itIntrinsic = sfm_data.intrinsics.begin();
        itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
      {
        const IndexT indexCam = itIntrinsic->first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        itIntrinsic->second.get()->updateFromParams(vec_params);
      }
    }
    // Structure is already updated directly if needed (no data wrapping)
    return true;
  }
}

bool Bundle_Adjustment_Ceres::EstimateUncertainty
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const Optimize_Options options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
 for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
  {
    const IndexT indexPose = itPose->first;

    const Pose3 & pose = itPose->second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
    itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
  {
    const IndexT indexCam = itIntrinsic->first;

    if (isValid(itIntrinsic->second->getType()))
    {
      map_intrinsics[indexCam] = itIntrinsic->second->getParams();

      double * parameter_block = &map_intrinsics[indexCam][0];
      problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
      if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
      {
        // set the whole parameter block as constant for best performance
        problem.SetParameterBlockConstant(parameter_block);
      }
      else
      {
        const std::vector<int> vec_constant_intrinsic =
          itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
        if (!vec_constant_intrinsic.empty())
        {
          ceres::SubsetParameterization *subset_parameterization =
            new ceres::SubsetParameterization(
              map_intrinsics[indexCam].size(), vec_constant_intrinsic);
          problem.SetParameterization(parameter_block, subset_parameterization);
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options_.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks!= sfm_data.structure.end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;

    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(itObs->first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

      if (cost_function)
        problem.AddResidualBlock(cost_function,
          p_LossFunction,
          &map_intrinsics[view->id_intrinsic][0],
          &map_poses[view->id_pose][0],
          iterTracks->second.X.data());
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(iterTracks->second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (Landmarks::iterator iterGCPTracks = sfm_data.control_points.begin();
      iterGCPTracks!= sfm_data.control_points.end(); ++iterGCPTracks)
    {
      const Observations & obs = iterGCPTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            itObs->second.x,
            options.control_point_opt.weight);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            nullptr,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            iterGCPTracks->second.X.data());
      }
      // Set the 3D point as FIXED (it's a GCP)
      problem.SetParameterBlockConstant(iterGCPTracks->second.X.data());
    }
  }

  // Basic sizes of parameters
  const int n_cam_ext_param = 6;

  // Confgure problem evaluator
  ceres::Problem::EvaluateOptions evaluate_options_;
  evaluate_options_.num_threads = ceres_options_.nb_threads_;
  ceres::CRSMatrix jacobian;
  double cost;
  // Evaluate problem
  problem.Evaluate(evaluate_options_, &cost, NULL, NULL, &jacobian);

  const int num_J_rows = jacobian.num_rows;
  const int num_J_cols = jacobian.num_cols;
  const int num_J_nonzeros = jacobian.values.size();

  const int num_cam_ext_param = n_cam_ext_param*sfm_data.poses.size();
  const int num_landmark_param = 3*sfm_data.structure.size();
  const int num_control_param = 3*sfm_data.control_points.size();

  // Determine the number of total intrinsic parameters
  int num_intrinsic_param = 0;
  Hash_Map<IndexT, int> start_intrinsic_param_per_row;
  int max_intrinsic_param_per_row=0;
  for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
    itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
  {
    const int n_intrinsic_param = map_intrinsics[itIntrinsic->first].size();
    // Save the start index of the intrinsics
    start_intrinsic_param_per_row[itIntrinsic->first] = num_intrinsic_param;
    // Compute total number of intrinsics
    num_intrinsic_param += n_intrinsic_param;
    // Compute max number of intrinsics
    if(n_intrinsic_param>max_intrinsic_param_per_row)
      max_intrinsic_param_per_row = n_intrinsic_param;
  }

  const int num_A_B_rows = num_J_rows;
  const int num_A_cols = num_cam_ext_param + num_intrinsic_param;
  const int num_B_cols = num_landmark_param + num_control_param;

  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> EigenSparseMatrix;
  // Convert the matrix to column major order
  EigenSparseMatrix sparse_jacobian =
      Eigen::MappedSparseMatrix<double, Eigen::RowMajor>(
          jacobian.num_rows, jacobian.num_cols,
          static_cast<int>(jacobian.values.size()),
          jacobian.rows.data(), jacobian.cols.data(), jacobian.values.data());

  // Matrix A (camera part of Jacobian) - each row has only values for one camera+intrinsics
  EigenSparseMatrix sparse_J_A = sparse_jacobian.block(0,0,num_A_B_rows,num_cam_ext_param+num_intrinsic_param);

  EigenSparseMatrix ExAll(num_A_B_rows,num_A_B_rows);
  for (int i=0;i<num_A_B_rows;i++){
    ExAll.insert(i,i) = 0.5;
  }

  EigenSparseMatrix SS = sparse_J_A.transpose()*ExAll*sparse_J_A;


  // Matrix B (landmark part of Jacobian) - each row has only values for one point
  EigenSparseMatrix sparse_J_B = sparse_jacobian.block(0,num_A_cols,num_A_B_rows,num_B_cols);

  // Find appropriate observations for each camera
  typedef Hash_Map< IndexT, Observations > ObservationsPerCam;
  ObservationsPerCam observations_per_cam;
  size_t global_obs_id = 0;

  for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks!= sfm_data.structure.end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {

      (observations_per_cam[itObs->first])[global_obs_id] = itObs->second;
      global_obs_id++;

    }
  }
  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (Landmarks::iterator iterGCPTracks = sfm_data.control_points.begin();
      iterGCPTracks!= sfm_data.control_points.end(); ++iterGCPTracks)
    {
      const Observations & obs = iterGCPTracks->second.obs;
      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        (observations_per_cam[itObs->first])[global_obs_id] = itObs->second;
        global_obs_id++;
      }
    }
  }

  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(num_cam_ext_param + num_intrinsic_param,num_cam_ext_param + num_intrinsic_param);
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(num_landmark_param + num_control_param,num_landmark_param + num_control_param);
  Eigen::MatrixXd camBlockMatrix = Eigen::MatrixXd::Zero(n_cam_ext_param,n_cam_ext_param);


  Eigen::Matrix2d Ex;
  Ex << 0.5,0,0,0.5;
  Eigen::MatrixXd J_cam_matrix;
  Eigen::MatrixXd J_intrinsics_matrix;
  // A' * Ex * A (cams)
  for (ObservationsPerCam::iterator iterObsCam = observations_per_cam.begin();
      iterObsCam!= observations_per_cam.end(); ++iterObsCam)
  {

    const IndexT view_id = iterObsCam->first;
    const IndexT pose_id = sfm_data.views.at(view_id)->id_pose;
    const IndexT intrinsic_id = sfm_data.views.at(view_id)->id_intrinsic;
    const int n_intrinsic_param = sfm_data.intrinsics.at(intrinsic_id)->getParams().size();

    camBlockMatrix.setZero();
    Eigen::MatrixXd intrinsicsCamBlockMatrix = Eigen::MatrixXd::Zero(n_cam_ext_param,n_intrinsic_param);
    Eigen::MatrixXd intrinsicsBlockMatrix = Eigen::MatrixXd::Zero(n_intrinsic_param,n_intrinsic_param);

    const Observations & obs = iterObsCam->second;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const IndexT obsID = itObs->first;
      // Get block from A
      J_cam_matrix = sparse_jacobian.block(obsID*2,pose_id*n_cam_ext_param,2,n_cam_ext_param);
      J_intrinsics_matrix = sparse_jacobian.block(obsID*2,num_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],2,n_intrinsic_param);

      camBlockMatrix += J_cam_matrix.transpose()*Ex*J_cam_matrix;
      intrinsicsCamBlockMatrix += J_cam_matrix.transpose()*Ex*J_intrinsics_matrix;
      intrinsicsBlockMatrix += J_intrinsics_matrix.transpose()*Ex*J_intrinsics_matrix;
    }
    // Pose diagonals
    U.block(pose_id*n_cam_ext_param,pose_id*n_cam_ext_param,n_cam_ext_param,n_cam_ext_param) = camBlockMatrix;
    // Intrinsic diagonals
    U.block(num_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],num_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],n_intrinsic_param,n_intrinsic_param) += intrinsicsBlockMatrix;
    // Pose - intrinsic diagonals
    U.block(pose_id*n_cam_ext_param,num_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],n_cam_ext_param,n_intrinsic_param) = intrinsicsCamBlockMatrix;
    U.block(num_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],pose_id*n_cam_ext_param,n_intrinsic_param,n_cam_ext_param) = intrinsicsCamBlockMatrix.transpose();
  }







  if (ceres_options_.bVerbose_)
  {
    // Display statistics about the minimization
    std::cout << std::endl
      << "After Bundle Adjustment estimation statistics (approximated RMSE):\n"
      << " #views: " << sfm_data.views.size() << "\n"
      << " #poses: " << sfm_data.poses.size() << "\n"
      << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
      << " #tracks: " << sfm_data.structure.size() << "\n"
      << " #J rows: " << num_J_rows << "\n"
      << " #J cols: " << num_J_cols << "\n"
      << " #J NonZero: " << global_obs_id << "\n"
      << " #Sparse J rows: " << sparse_jacobian.rows()<<" :: "<<sparse_jacobian.cols() << "\n"
      << " #Sparse JA rows: " << sparse_J_A.rows()<<" :: "<<sparse_J_A.cols() << "\n"
      << " #Sparse JB rows: " << sparse_J_B.rows()<<" :: "<<sparse_J_B.cols() << "\n"
      << " #Sparse J nonZero: " << sparse_jacobian.nonZeros() << "\n"
      << " #Sparse JA nonZero: " << sparse_J_A.nonZeros() << "\n"
      << " #Sparse JB nonZero: " << sparse_J_B.nonZeros() << "\n"
      << " Final RMSE: " << std::sqrt( cost / num_J_rows) << "\n"
      << " Final cost: " << cost << "\n"
      << " --------------------------\n"
      << " U: \n" << U << "\n"
      << " UU: \n "<<SS<<"\n"
      << std::endl;
  }


  return true;
}


} // namespace sfm
} // namespace openMVG

