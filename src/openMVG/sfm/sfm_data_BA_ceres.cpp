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
  const int single_cam_ext_param = 6;

  // Confgure problem evaluator
  ceres::Problem::EvaluateOptions evaluate_options_;
  evaluate_options_.num_threads = ceres_options_.nb_threads_;
  ceres::CRSMatrix jacobian;
  double cost;
  // Evaluate problem
  problem.Evaluate(evaluate_options_, &cost, NULL, NULL, &jacobian);

  // -----------------------------------------------
  // Compute size of parameters
  // -----------------------------------------------
  
  // Compute the size of Jacobian matrix
  const int num_J_rows = jacobian.num_rows;
  const int num_J_cols = jacobian.num_cols;
  const int num_J_nonzeros = jacobian.values.size();
  // Compute the size of different parameter blocks
  const int total_cam_ext_param = single_cam_ext_param*sfm_data.poses.size();
  const int total_landmark_param = 3*sfm_data.structure.size();
  const int total_control_param = 3*sfm_data.control_points.size();
  int total_intrinsic_param = 0;
  int max_intrinsic_param = 0;
  // Save the index of the start of next intrinsic block of params
  Hash_Map<IndexT, int> start_intrinsic_param_per_row;
  for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
    itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
  {
    const int n_intrinsic_param = map_intrinsics[itIntrinsic->first].size();
    // Save the start index of the intrinsics
    start_intrinsic_param_per_row[itIntrinsic->first] = total_intrinsic_param;
    // Compute total number of intrinsics
    total_intrinsic_param += n_intrinsic_param;
    // Compute max number of intrinsics
    if(n_intrinsic_param>max_intrinsic_param)
      max_intrinsic_param = n_intrinsic_param;
  }

  // -----------------------------------------------
  // Convert Ceres Jacobian to Sparse Eigen Matrix
  // -----------------------------------------------
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor> EigenSparseMatrix;
  // Convert the matrix to column major order
  EigenSparseMatrix sparse_jacobian =
      Eigen::MappedSparseMatrix<double, Eigen::RowMajor>(
          jacobian.num_rows, jacobian.num_cols,
          static_cast<int>(jacobian.values.size()),
          jacobian.rows.data(), jacobian.cols.data(), jacobian.values.data());


  // -----------------------------------------------
  // Determine problem dependent IDs
  // -----------------------------------------------
  typedef Hash_Map< IndexT, std::shared_ptr<Observation> > Problem_ObservationID;
  Problem_ObservationID problem_observationID;
  typedef Hash_Map< IndexT, std::shared_ptr<Landmark> > Problem_LandmarkID;
  Problem_LandmarkID problem_landmarkID;
  typedef Hash_Map< IndexT, std::vector<IndexT> > Problem_ObservationsPerView;
  Problem_ObservationsPerView problem_obs_per_view;
  typedef Hash_Map< IndexT, Hash_Map<IndexT, IndexT> > Problem_ObservationPerLandmarkView;
  Problem_ObservationPerLandmarkView problem_obs_per_land_view;
  
  Eigen::VectorXi nonZeroW = Eigen::VectorXi::Zero(total_cam_ext_param + total_intrinsic_param);

  IndexT pObsID=0;
  IndexT pLandID=0;
  for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks!= sfm_data.structure.end(); ++iterTracks)
  {
    std::set<IndexT> addedIntrinsics;
    const Observations & obs = iterTracks->second.obs;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      problem_observationID[pObsID] = std::make_shared<Observation>(itObs->second);
      (problem_obs_per_view[itObs->first]).push_back(pObsID);
      (problem_obs_per_land_view[pLandID])[itObs->first] = pObsID;

	  const IndexT view_id = itObs->first;
	  const IndexT pose_id = sfm_data.views.at(view_id)->id_pose;
      const IndexT intrinsic_id = sfm_data.views.at(view_id)->id_intrinsic;
      const int n_intrinsic_param = sfm_data.intrinsics.at(intrinsic_id)->getParams().size();
      
      for (int i=0;i<3;i++){
		nonZeroW(pose_id*3+i) += 3;		 
      }
      
      if (addedIntrinsics.find(intrinsic_id) == addedIntrinsics.end())
      {
        // Point was not added for this intrinsics group yet
        for (int i=0;i<n_intrinsic_param;i++){
		  nonZeroW(total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id] + i) += 3;
        }
        addedIntrinsics.insert(intrinsic_id);
      }

      pObsID++;
    }
    problem_landmarkID[pLandID] = std::make_shared<Landmark>(iterTracks->second);
    pLandID++;
  }
  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (Landmarks::iterator iterGCPTracks = sfm_data.control_points.begin();
      iterGCPTracks!= sfm_data.control_points.end(); ++iterGCPTracks)
    {
      std::set<IndexT> addedIntrinsics;
      const Observations & obs = iterGCPTracks->second.obs;
      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        problem_observationID[pObsID] = std::make_shared<Observation>(itObs->second);
        (problem_obs_per_view[itObs->first]).push_back(pObsID);
        (problem_obs_per_land_view[pLandID])[itObs->first] = pObsID;
        
        const IndexT view_id = itObs->first;
	    const IndexT pose_id = sfm_data.views.at(view_id)->id_pose;
        const IndexT intrinsic_id = sfm_data.views.at(view_id)->id_intrinsic;
        const int n_intrinsic_param = sfm_data.intrinsics.at(intrinsic_id)->getParams().size();
      
        for (int i=0;i<3;i++){
		  nonZeroW(pose_id*3+i) += 3;		 
        }
      
        if (addedIntrinsics.find(intrinsic_id) == addedIntrinsics.end())
        {
          // Point was not added for this intrinsics group yet
          for (int i=0;i<n_intrinsic_param;i++){
		    nonZeroW(total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id] + i) += 3;
          }
          addedIntrinsics.insert(intrinsic_id);
        }
        pObsID++;
      }
      problem_landmarkID[pLandID] = std::make_shared<Landmark>(iterGCPTracks->second);
      pLandID++;  
    }  
  }

  // -----------------------------------------------
  // Result matrices
  // -----------------------------------------------
  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(total_cam_ext_param + total_intrinsic_param,total_cam_ext_param + total_intrinsic_param);
  EigenSparseMatrix V_inv(total_landmark_param + total_control_param,total_landmark_param + total_control_param);
  EigenSparseMatrix W(total_cam_ext_param + total_intrinsic_param,total_landmark_param + total_control_param);
  Eigen::MatrixXd WUW = Eigen::MatrixXd::Zero(total_cam_ext_param + total_intrinsic_param,total_cam_ext_param + total_intrinsic_param);
  //EigenSparseMatrix WVW(total_cam_ext_param + total_intrinsic_param,total_cam_ext_param + total_intrinsic_param);

  // Uncertainty of each feature detected
  Eigen::Matrix2d Ex;
  Ex << 0.5,0,0,0.5;
  
  
  // -----------------------------------------------
  // Compute U
  // -----------------------------------------------
  Eigen::MatrixXd camBlockMatrix = Eigen::MatrixXd::Zero(single_cam_ext_param,single_cam_ext_param);
  Eigen::MatrixXd J_cam_matrix;
  Eigen::MatrixXd J_intrinsics_matrix;
  
  for (Problem_ObservationsPerView::iterator iterObsView = problem_obs_per_view.begin();
      iterObsView!= problem_obs_per_view.end(); ++iterObsView)
  {
    const IndexT view_id = iterObsView->first;
    const IndexT pose_id = sfm_data.views.at(view_id)->id_pose;
    const IndexT intrinsic_id = sfm_data.views.at(view_id)->id_intrinsic;
    const int n_intrinsic_param  = sfm_data.intrinsics.at(intrinsic_id)->getParams().size();
    // Initialize sum computations
    camBlockMatrix.setZero();
    Eigen::MatrixXd intrinsicsCamBlockMatrix = Eigen::MatrixXd::Zero(single_cam_ext_param,n_intrinsic_param);
    Eigen::MatrixXd intrinsicsBlockMatrix = Eigen::MatrixXd::Zero(n_intrinsic_param,n_intrinsic_param);
    
    // Loop through observations and comptue sum of pairwise blocks
    const std::vector<IndexT> obs_ids = iterObsView->second;
    for (std::vector<IndexT>::const_iterator itObs = obs_ids.begin();
      itObs != obs_ids.end(); ++itObs)
    {
	  const IndexT obsID = *itObs;
	  // Get blocks from Jacobian matrix
	  J_cam_matrix = sparse_jacobian.block(obsID*2,pose_id*single_cam_ext_param,2,single_cam_ext_param);
      J_intrinsics_matrix = sparse_jacobian.block(obsID*2,total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],2,n_intrinsic_param);
		
	  // Compute resulting blocks
      camBlockMatrix += J_cam_matrix.transpose()*Ex*J_cam_matrix;
      intrinsicsCamBlockMatrix += J_cam_matrix.transpose()*Ex*J_intrinsics_matrix;
      intrinsicsBlockMatrix += J_intrinsics_matrix.transpose()*Ex*J_intrinsics_matrix;
    }
    // Save results
    // Pose diagonals
    U.block(pose_id*single_cam_ext_param,pose_id*single_cam_ext_param,single_cam_ext_param,single_cam_ext_param) = camBlockMatrix;
    // Intrinsic diagonals
    U.block(total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],n_intrinsic_param,n_intrinsic_param) += intrinsicsBlockMatrix;
    // Pose - intrinsic diagonals
    U.block(pose_id*single_cam_ext_param,total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],single_cam_ext_param,n_intrinsic_param) = intrinsicsCamBlockMatrix;
    U.block(total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],pose_id*single_cam_ext_param,n_intrinsic_param,single_cam_ext_param) = intrinsicsCamBlockMatrix.transpose();    
  }
  
  
  // -----------------------------------------------
  // Compute V_inverse and W
  // -----------------------------------------------
  Eigen::MatrixXd pointBlockMatrix = Eigen::MatrixXd::Zero(3,3);
  Eigen::MatrixXd camPointBlockMatrix = Eigen::MatrixXd::Zero(single_cam_ext_param,3);
  Eigen::MatrixXd J_point_matrix;
  // Reserve space in V: Each row will have 3 non-zero elements
  V_inv.reserve(3);
  W.reserve(nonZeroW);

  // Loop through landmarks
  for (Problem_ObservationPerLandmarkView::iterator iterTracks = problem_obs_per_land_view.begin();
    iterTracks!= problem_obs_per_land_view.end(); ++iterTracks)
  {
    const IndexT track_id = iterTracks->first;
    // Initialize computation
    pointBlockMatrix.setZero();
    const Hash_Map<IndexT,IndexT> * obs_per_land_view = &(iterTracks->second);
    for (Hash_Map<IndexT,IndexT>::const_iterator iterObsLV = obs_per_land_view->begin();
	  iterObsLV!=obs_per_land_view->end();++iterObsLV){
	  const IndexT view_id = iterObsLV->first;
	  const IndexT pose_id = sfm_data.views.at(view_id)->id_pose;
      const IndexT intrinsic_id = sfm_data.views.at(view_id)->id_intrinsic;
      const int n_intrinsic_param = sfm_data.intrinsics.at(intrinsic_id)->getParams().size();
	  const IndexT obs_id = iterObsLV->second;
	  
	  // Get track/obs block from Jacobian matrix
	  J_cam_matrix = sparse_jacobian.block(obs_id*2,pose_id*single_cam_ext_param,2,single_cam_ext_param);
      J_intrinsics_matrix = sparse_jacobian.block(obs_id*2,total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],2,n_intrinsic_param);
      J_point_matrix = sparse_jacobian.block(obs_id*2,(total_cam_ext_param + total_intrinsic_param)+track_id*3,2,3);
      
      
      // Compute diagonal block for V_inverse
      pointBlockMatrix += J_point_matrix.transpose()*Ex*J_point_matrix;
	  
	  // Compute (cam,landmark) block for W
	  camPointBlockMatrix = J_cam_matrix.transpose()*Ex*J_point_matrix;
	  // Compute (intrinsic,landmark) block for W
	  Eigen::MatrixXd intrinsicsPointBlockMatrix = J_intrinsics_matrix.transpose()*Ex*J_point_matrix;
	  
	  // Save block to W
	  for(int r=0;r<single_cam_ext_param;r++){
        for(int c=0;c<3;c++){
          W.insert(pose_id*single_cam_ext_param+r,track_id*3+c) = camPointBlockMatrix(r,c);
        }
      }
      
	  for(int r=0;r<n_intrinsic_param;r++){
        for(int c=0;c<3;c++){
          W.coeffRef(total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id]+r,track_id*3+c) += intrinsicsPointBlockMatrix(r,c);
        }
      }
	  
	}
	
    // Compute inverse of block of V
    pointBlockMatrix=pointBlockMatrix.inverse();
    
    // Save block to V
	for(int r=0;r<3;r++){
      for(int c=0;c<3;c++){
        V_inv.insert(track_id*3+r,track_id*3+c) = pointBlockMatrix(r,c);
      }
    }    
  }
  
  // -----------------------------------------------
  // WVW
  // -----------------------------------------------
  Eigen::MatrixXd camWBlockMatrix;
  Eigen::MatrixXd camWVWBlockMatrix;
  // Loop through landmarks
  for (Problem_ObservationPerLandmarkView::iterator iterTracks = problem_obs_per_land_view.begin();
    iterTracks!= problem_obs_per_land_view.end(); ++iterTracks)
  {
    const IndexT track_id = iterTracks->first;
    // Initialize computation
    pointBlockMatrix.setZero();
    const Hash_Map<IndexT,IndexT> * obs_per_land_view = &(iterTracks->second);
    for (Hash_Map<IndexT,IndexT>::const_iterator iterObsLV = obs_per_land_view->begin();
	  iterObsLV!=obs_per_land_view->end();++iterObsLV){
	  const IndexT view_id = iterObsLV->first;
	  const IndexT pose_id = sfm_data.views.at(view_id)->id_pose;
      const IndexT intrinsic_id = sfm_data.views.at(view_id)->id_intrinsic;
      const int n_intrinsic_param = sfm_data.intrinsics.at(intrinsic_id)->getParams().size();
	  const IndexT obs_id = iterObsLV->second;
	  
	  // Get track/obs block from Jacobian matrix
	  camWBlockMatrix = W.block(pose_id*single_cam_ext_param,track_id*3,single_cam_ext_param,3);
	  camWVWBlockMatrix = camWBlockMatrix*V_inv.block(track_id*3,track_id*3,3,3)*camWBlockMatrix.transpose();
	  
	  WUW.block(track_id*3,track_id*3,3,3) = camWBlockMatrix;
	  
	}	   
  }
  
  Eigen::MatrixXd UWVW = ((U-WUW).transpose()).inverse();
  
  

/*
  int gObsID=0;
  int gTrackID=0;
  // Loop through landmarks
  for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks!= sfm_data.structure.end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    // Set the resulting block to zero
    pointBlockMatrix.setZero();

    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const IndexT view_id = itObs->first;
      const IndexT pose_id = sfm_data.views.at(view_id)->id_pose;
      const IndexT intrinsic_id = sfm_data.views.at(view_id)->id_intrinsic;
      const int n_intrinsic_param = sfm_data.intrinsics.at(intrinsic_id)->getParams().size();

      // Get B matrix from J
      J_point_matrix = sparse_jacobian.block(gObsID*2,(total_cam_ext_param + total_intrinsic_param)+gTrackID*3,2,3);
      pointBlockMatrix += J_point_matrix.transpose()*Ex*J_point_matrix;

      J_A_matrix = sparse_jacobian.block(gObsID*2,pose_id*single_cam_ext_param,2,single_cam_ext_param);
      J_B_matrix = sparse_jacobian.block(gObsID*2,(total_cam_ext_param + total_intrinsic_param)+gTrackID*3,2,3);
      camPointBlockMatrix = J_A_matrix.transpose()*Ex*J_B_matrix;

      for(int r=0;r<3;r++){
        for(int c=0;c<3;c++){
          W.insert(pose_id*single_cam_ext_param*3+r,gTrackID*3+c) = camPointBlockMatrix(r,c);
        }
      }



      gObsID++;
    }
    pointBlockMatrix=pointBlockMatrix.inverse();

    for(int r=0;r<3;r++){
      for(int c=0;c<3;c++){
        V_inv.insert(gTrackID*3+r,gTrackID*3+c) = pointBlockMatrix(r,c);
      }
    }
    gTrackID++;
  }
  
  */

  // -----------------------------------------------
  // Compute Partial Jacobian matrices (Eigen)
  // -----------------------------------------------
  const int num_A_B_rows = num_J_rows;
  const int num_A_cols = total_cam_ext_param + total_intrinsic_param;
  const int num_B_cols = total_landmark_param + total_control_param;
  // Matrix A (camera part of Jacobian) - each row has only values for one camera+intrinsics
  EigenSparseMatrix sparse_J_A = sparse_jacobian.block(0,0,num_A_B_rows,total_cam_ext_param+total_intrinsic_param);
  // Matrix B (landmark part of Jacobian) - each row has only values for one point
  EigenSparseMatrix sparse_J_B = sparse_jacobian.block(0,num_A_cols,num_A_B_rows,num_B_cols);

  EigenSparseMatrix ExAll(num_A_B_rows,num_A_B_rows);
  for (int i=0;i<num_A_B_rows;i++){
    ExAll.insert(i,i) = 0.5;
  }
  
  EigenSparseMatrix JJA = sparse_J_A.transpose()*ExAll*sparse_J_A;
  EigenSparseMatrix JJB = sparse_J_B.transpose()*ExAll*sparse_J_B;
  EigenSparseMatrix JJW = sparse_J_A.transpose()*ExAll*sparse_J_B;
  
  EigenSparseMatrix JJWVW = JJW*JJB*JJW.transpose();
  Eigen::MatrixXd UJJWUW = (U-Eigen::MatrixXd(JJWVW));
  UJJWUW = UJJWUW.transpose().inverse();

/*
  // Observations for each Landmark and camera
  typedef Hash_Map< IndexT, Hash_Map< IndexT, Observations > > ObservationsPerLandmarkAndCam;
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

  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(total_cam_ext_param + total_intrinsic_param,total_cam_ext_param + total_intrinsic_param);
  EigenSparseMatrix V_inv(total_landmark_param + total_control_param,total_landmark_param + total_control_param);
  EigenSparseMatrix W(total_cam_ext_param + total_intrinsic_param,total_landmark_param + total_control_param);


  Eigen::MatrixXd camBlockMatrix = Eigen::MatrixXd::Zero(single_cam_ext_param,single_cam_ext_param);

  // Uncertainty of each feature detected
  Eigen::Matrix2d Ex;
  Ex << 0.5,0,0,0.5;

  Eigen::MatrixXd J_cam_matrix;
  Eigen::MatrixXd J_point_matrix;
  Eigen::MatrixXd J_intrinsics_matrix;
  Eigen::MatrixXd J_A_matrix;
  Eigen::MatrixXd J_B_matrix;
  // A' * Ex * A (cams)
  for (ObservationsPerCam::iterator iterObsCam = observations_per_cam.begin();
      iterObsCam!= observations_per_cam.end(); ++iterObsCam)
  {

    const IndexT view_id = iterObsCam->first;
    const IndexT pose_id = sfm_data.views.at(view_id)->id_pose;
    const IndexT intrinsic_id = sfm_data.views.at(view_id)->id_intrinsic;
    const int n_intrinsic_param = sfm_data.intrinsics.at(intrinsic_id)->getParams().size();

    camBlockMatrix.setZero();
    Eigen::MatrixXd intrinsicsCamBlockMatrix = Eigen::MatrixXd::Zero(single_cam_ext_param,n_intrinsic_param);
    Eigen::MatrixXd intrinsicsBlockMatrix = Eigen::MatrixXd::Zero(n_intrinsic_param,n_intrinsic_param);

    const Observations & obs = iterObsCam->second;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const IndexT obsID = itObs->first;
      // Get block from A
      J_cam_matrix = sparse_jacobian.block(obsID*2,pose_id*single_cam_ext_param,2,single_cam_ext_param);
      J_intrinsics_matrix = sparse_jacobian.block(obsID*2,total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],2,n_intrinsic_param);

      camBlockMatrix += J_cam_matrix.transpose()*Ex*J_cam_matrix;
      intrinsicsCamBlockMatrix += J_cam_matrix.transpose()*Ex*J_intrinsics_matrix;
      intrinsicsBlockMatrix += J_intrinsics_matrix.transpose()*Ex*J_intrinsics_matrix;
    }
    // Pose diagonals
    U.block(pose_id*single_cam_ext_param,pose_id*single_cam_ext_param,single_cam_ext_param,single_cam_ext_param) = camBlockMatrix;
    // Intrinsic diagonals
    U.block(total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],n_intrinsic_param,n_intrinsic_param) += intrinsicsBlockMatrix;
    // Pose - intrinsic diagonals
    U.block(pose_id*single_cam_ext_param,total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],single_cam_ext_param,n_intrinsic_param) = intrinsicsCamBlockMatrix;
    U.block(total_cam_ext_param + start_intrinsic_param_per_row[intrinsic_id],pose_id*single_cam_ext_param,n_intrinsic_param,single_cam_ext_param) = intrinsicsCamBlockMatrix.transpose();
  }

  Eigen::MatrixXd pointBlockMatrix = Eigen::MatrixXd::Zero(3,3);
  Eigen::MatrixXd camPointBlockMatrix = Eigen::MatrixXd::Zero(single_cam_ext_param,3);
  // Each row in V will have 3 non-zero elements
  V_inv.reserve(3);


  int gObsID=0;
  int gTrackID=0;
  // Loop through landmarks
  for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks!= sfm_data.structure.end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    // Set the resulting block to zero
    pointBlockMatrix.setZero();

    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const IndexT view_id = itObs->first;
      const IndexT pose_id = sfm_data.views.at(view_id)->id_pose;
      const IndexT intrinsic_id = sfm_data.views.at(view_id)->id_intrinsic;
      const int n_intrinsic_param = sfm_data.intrinsics.at(intrinsic_id)->getParams().size();

      // Get B matrix from J
      J_point_matrix = sparse_jacobian.block(gObsID*2,(total_cam_ext_param + total_intrinsic_param)+gTrackID*3,2,3);
      pointBlockMatrix += J_point_matrix.transpose()*Ex*J_point_matrix;

      J_A_matrix = sparse_jacobian.block(gObsID*2,pose_id*single_cam_ext_param,2,single_cam_ext_param);
      J_B_matrix = sparse_jacobian.block(gObsID*2,(total_cam_ext_param + total_intrinsic_param)+gTrackID*3,2,3);
      camPointBlockMatrix = J_A_matrix.transpose()*Ex*J_B_matrix;

      for(int r=0;r<3;r++){
        for(int c=0;c<3;c++){
          W.insert(pose_id*single_cam_ext_param*3+r,gTrackID*3+c) = camPointBlockMatrix(r,c);
        }
      }



      gObsID++;
    }
    pointBlockMatrix=pointBlockMatrix.inverse();

    for(int r=0;r<3;r++){
      for(int c=0;c<3;c++){
        V_inv.insert(gTrackID*3+r,gTrackID*3+c) = pointBlockMatrix(r,c);
      }
    }
    gTrackID++;
  }




*/



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
      << " #Sparse J rows: " << sparse_jacobian.rows()<<" :: "<<sparse_jacobian.cols() << "\n"
      /*<< " #Sparse JA rows: " << sparse_J_A.rows()<<" :: "<<sparse_J_A.cols() << "\n"
      << " #Sparse JB rows: " << sparse_J_B.rows()<<" :: "<<sparse_J_B.cols() << "\n"
      << " #Sparse J nonZero: " << sparse_jacobian.nonZeros() << "\n"
      << " #Sparse JA nonZero: " << sparse_J_A.nonZeros() << "\n"
      << " #Sparse JB nonZero: " << sparse_J_B.nonZeros() << "\n"*/
      << " Final RMSE: " << std::sqrt( cost / num_J_rows) << "\n"
      << " Final cost: " << cost << "\n"
      << " --------------------------\n"
      /*<< " U: \n" << U << "\n"
      << " UU: \n "<<JJA<<"\n"*/
      /*<< " V: \n" << V_inv.block(100,100,10,10) << "\n"
      << " VV: \n "<<JJB.block(100,100,10,10)<<"\n"*/
      << " W: \n" << UWVW.block(0,0,W.rows(),10) << "\n"
      << " WW: \n "<<UJJWUW.block(0,0,JJW.rows(),10)<<"\n"
      
      << " #Sparse W rows: " << W.rows()<<" :: "<<W.cols() << "\n"
      << " #Sparse JW rows: " << JJW.rows()<<" :: "<<JJW.cols() << "\n"
      /*<< " V: \n" << V_inv.block(100,100,10,10) << "\n"
      << " VV: \n "<<SS.block(100,100,10,10)<<"\n"*/
      << std::endl;
  }


  return true;
}


} // namespace sfm
} // namespace openMVG

