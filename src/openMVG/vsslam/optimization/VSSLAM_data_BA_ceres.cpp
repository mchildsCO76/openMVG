
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



#include <openMVG/vsslam/optimization/VSSLAM_data_BA_ceres.hpp>
#include <openMVG/cameras/Camera_Common.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA_ceres_camera_functor.hpp>
#include <openMVG/vsslam/mapping/VSSLAM_data_io_ply.hpp>

#include "openMVG/geometry/Similarity3_Kernel.hpp"
//- Robust estimation - LMeds (since no threshold can be defined)
#include "openMVG/robust_estimation/robust_estimator_LMeds.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/types.hpp"

using namespace openMVG;

namespace openMVG {
namespace VSSLAM {

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
      return ResidualErrorFunctor_Pinhole_Intrinsic_Rts::Create(observation, weight);
     break;
    default:
      return nullptr;
  }
}


VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options::BA_Ceres_options
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

VSSLAM_Bundle_Adjustment_Ceres::VSSLAM_Bundle_Adjustment_Ceres
(
  VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options options
)
: ceres_options_(options)
{}

VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options &
VSSLAM_Bundle_Adjustment_Ceres::ceres_options()
{
  return ceres_options_;
}


bool VSSLAM_Bundle_Adjustment_Ceres::Adjust
(
  VSSLAM_Data & slam_data,     // the SfM scene to refine
  const sfm::Optimize_Options options
)
{
  //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------

  Save_PLY(slam_data,"before.ply",sfm::ESfM_Data::ALL);
  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & kF_it : slam_data.keyframes)
  {
    const size_t kF_idx = kF_it.first;

    const Similarity3 & pose = kF_it.second->pose_;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();
    const double s = pose.scale();


    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[kF_idx] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2),s};

    double * parameter_block = &map_poses[kF_idx][0];
    problem.AddParameterBlock(parameter_block, 7);
    if (options.extrinsics_opt == sfm::Extrinsic_Parameter_Type::NONE)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == sfm::Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == sfm::Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
        vec_constant_extrinsic.push_back(6);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(7, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : slam_data.cam_intrinsics)
  {
    const IndexT cam_idx = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[cam_idx] = intrinsic_it.second->getParams();

      double * parameter_block = &map_intrinsics[cam_idx][0];
      problem.AddParameterBlock(parameter_block, map_intrinsics[cam_idx].size());
      if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
      {
        // set the whole parameter block as constant for best performance
        problem.SetParameterBlockConstant(parameter_block);
      }
      else
      {
        const std::vector<int> vec_constant_intrinsic =
          intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
        if (!vec_constant_intrinsic.empty())
        {
          ceres::SubsetParameterization *subset_parameterization =
            new ceres::SubsetParameterization(
              map_intrinsics[cam_idx].size(), vec_constant_intrinsic);
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
  for (auto & structure_landmark_it : slam_data.structure)
  {
    const LandmarkObservations & obs = structure_landmark_it.second.obs_;

    for (const auto & obs_it : obs)
    {
      Frame * obs_frame = obs_it.first;
      // Build the residual block corresponding to the track observation:

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(slam_data.cam_intrinsics[obs_frame->camId_], obs_frame->getFeaturePosition(obs_it.second));
      if (cost_function)
        problem.AddResidualBlock(cost_function,
          p_LossFunction,
          &map_intrinsics[obs_frame->camId_][0],
          &map_poses[obs_frame->frameId_][0],
          structure_landmark_it.second.pt_.data());
    }
    if (options.structure_opt == sfm::Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.pt_.data());
  }


  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type = ceres_options_.preconditioner_type_;
  ceres_config_options.linear_solver_type = ceres_options_.linear_solver_type_;
  ceres_config_options.sparse_linear_algebra_library_type = ceres_options_.sparse_linear_algebra_library_type_;
  ceres_config_options.minimizer_progress_to_stdout = ceres_options_.bVerbose_;
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

    Save_PLY(slam_data,"after_nocam.ply",sfm::ESfM_Data::ALL);
    if (ceres_options_.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #keyframes: " << slam_data.keyframes.size() << "\n"
        << " #intrinsics: " << slam_data.cam_intrinsics.size() << "\n"
        << " #points: " << slam_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
    }
  }


  return false;
}

}
}
