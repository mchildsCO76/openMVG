
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



#include <openMVG/vsslam/optimization/VSSLAM_data_BA_ceres.hpp>
#include <openMVG/cameras/Camera_Common.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA_ceres_camera_functor.hpp>

#include "openMVG/geometry/Similarity3_Kernel.hpp"
//- Robust estimation - LMeds (since no threshold can be defined)
#include "openMVG/robust_estimation/robust_estimator_LMeds.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_transform.hpp"
#include "openMVG/types.hpp"

#include <openMVG/vsslam/tracking/PoseEstimation.hpp>
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
  VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options options = VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options()
)
:ceres_options_(options)
{}

VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options &
VSSLAM_Bundle_Adjustment_Ceres::ceres_options()
{
  return ceres_options_;
}

// Optimize free poses in the vector (not active -> not in global map)
// frame_i is the pose which is used for 3D-2D matching
// 3D Points in matches_3D_ptr_cur_idx remain freezed
// 3D Points in vec_triangulated_pts are optimized with all free poses
bool VSSLAM_Bundle_Adjustment_Ceres::OptimizePose
(
  std::vector<Frame*> * vec_frames,
  Frame * frame_i,
  Hash_Map<MapLandmark *,size_t> * matches_3D_ptr_cur_idx,
  std::vector<std::unique_ptr<MapLandmark> > * vec_triangulated_pts
)
{
  std::cout<<"Pose optimization: \n";
  // if vector of frames is empty and no frame_i is given -> no pose to optimize
  if ((!vec_frames || vec_frames->empty()) && !frame_i)
  {
    return false;
  }

  // If none points are
  /*if ((!matches_3D_ptr_cur_idx || matches_3D_ptr_cur_idx->empty()) && (!vec_triangulated_pts || vec_triangulated_pts->empty()))
  {
    return false;
  }*/

  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<size_t, std::vector<double> > map_intrinsics;
  Hash_Map<size_t, std::vector<double> > map_poses;

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));

  // Pose parameters
  Mat3 R;
  double angleAxis[3];
  Vec3 t;
  double s;


  std::cout<<"Pose optimization: Frames\n";
  // Add all poses in local scope
  // (Not active ones are not yet in global map and are free to move)
  if (vec_frames)
  {
    for (Frame * & frame : *vec_frames)
    {
      const size_t & frameId = frame->getFrameId();
      const size_t & cam_idx = frame->getCamId();

      // Get pose in the WORLD reference frame
      frame->getPose_cr_Rts(R,t,s);

      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      // angleAxis + translation
      map_poses[frameId] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2),s};
      double * parameter_block = &map_poses[frameId][0];
      problem.AddParameterBlock(parameter_block, 7);

      // Frame is active when its in the global map (cant be changed by localBA)
      if (frame->isActive())
        problem.SetParameterBlockConstant(parameter_block);

      if (map_intrinsics.find(cam_idx)==map_intrinsics.end())
      {
        map_intrinsics[cam_idx] = frame->getCameraIntrinsics()->getParams();
        double * parameter_block = &map_intrinsics[cam_idx][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[cam_idx].size());
        problem.SetParameterBlockConstant(parameter_block);
      }
    }
  }

  // Add frame of interest (current frame)
  // Matches 3D-2D are relating to this frame
  std::cout<<"Pose optimization: Frame\n";
  if (frame_i)
  {
    const size_t & frame_i_id = frame_i->getFrameId();
    const size_t & cam_i_idx = frame_i->getCamId();
    IntrinsicBase * & cam_i_intrinsic = frame_i->getCameraIntrinsics();

    // If pose is not already among the others added
    if (map_poses.find(frame_i_id) == map_poses.end())
    {
      // Get pose in the WORLD reference frame
      frame_i->getPose_cr_Rts(R,t,s);

      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      // angleAxis + translation
      map_poses[frame_i_id] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2),s};
      double * parameter_block = &map_poses[frame_i_id][0];
      problem.AddParameterBlock(parameter_block, 7);

      // Frame is active when its in the global map (cant be changed by localBA)
      if (frame_i->isActive())
        problem.SetParameterBlockConstant(parameter_block);

      if (map_intrinsics.find(cam_i_idx)==map_intrinsics.end())
      {
        map_intrinsics[cam_i_idx] = frame_i->getCameraIntrinsics()->getParams();
        double * parameter_block = &map_intrinsics[cam_i_idx][0];
        problem.AddParameterBlock(parameter_block, map_intrinsics[cam_i_idx].size());
        problem.SetParameterBlockConstant(parameter_block);
      }
    }

    // Impose 3D-2D matches (valid points are static) for frame_i
    if(matches_3D_ptr_cur_idx && !matches_3D_ptr_cur_idx->empty())
    {
      // Add all 2D-3D relations of the free pose
      for (auto & match : *matches_3D_ptr_cur_idx)
      {
        // Gets the position of the point (undistorted if distortion is known and calibrated or distorted if not calibrated)
        const Vec2 & obs_pt = frame_i->getFeaturePosition(match.second);

        // Build the residual block corresponding to the track observation:

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function = IntrinsicsToCostFunction(cam_i_intrinsic, obs_pt);

        if (cost_function)
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[cam_i_idx][0],
            &map_poses[frame_i_id][0],
            match.first->X_.data());
          // Set structure fixed if in the global map
          if (match.first->isActive())
            problem.SetParameterBlockConstant(match.first->X_.data());
        }
      }
    }

    // Impose 3D-2D matches that are already associated with frame_i (valid points are static)
    for (size_t mp_i = 0; mp_i < frame_i->map_points_.size(); ++mp_i)
    {
      MapLandmark* & map_point = frame_i->map_points_[mp_i];
      if (!map_point)
        continue;

      // Gets the position of the point (undistorted if distortion is known and calibrated or distorted if not calibrated)
      const Vec2 & obs_pt = frame_i->getFeaturePosition(mp_i);

      // Build the residual block corresponding to the track observation:

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction * cost_function = IntrinsicsToCostFunction(cam_i_intrinsic, obs_pt);


      if (cost_function)
      {
        problem.AddResidualBlock(cost_function,
          p_LossFunction,
          &map_intrinsics[cam_i_idx][0],
          &map_poses[frame_i_id][0],
          map_point->X_.data());
        // Set structure fixed if in the global map
        if (map_point->isActive())
          problem.SetParameterBlockConstant(map_point->X_.data());
      }
    }
  }

  std::cout<<"Pose optimization: Triangulated Pts\n";
  // Add newly triangulated points (3d points are free)
  if(vec_triangulated_pts && !vec_triangulated_pts->empty())
  {
    for (std::unique_ptr<MapLandmark> & point_it : *vec_triangulated_pts)
    {
      LandmarkObservations & obs = point_it->obs_;
      for(auto & m_o : obs)
      {
        Frame * frame_i = m_o.second.frame_ptr;
        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(frame_i->getCameraIntrinsics(), frame_i->getFeaturePosition(m_o.second.feat_id));

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[frame_i->getCamId()][0],
            &map_poses[frame_i->getFrameId()][0],
            point_it->X_.data());
      }
    }
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


  std::cout<<"Pose optimization: Solve\n";
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
        << " #keyframes: " << map_poses.size() << "\n"
        << " #intrinsics: " << map_intrinsics.size() << "\n"
        //<< " #points: " << slam_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
    }


    std::cout<<"Pose optimization: Updated poses\n";
    // Update camera poses with refined data
    Mat3 R_refined;
    Vec3 t_refined;
    double s_refined;
    if (vec_frames)
    {
      for (Frame * & frame : *vec_frames)
      {
        // Active frames are used in global map -> fixed
        if (frame->isActive())
          continue;
        // Get frame and update values
        size_t frame_id = frame->getFrameId();
        ceres::AngleAxisToRotationMatrix(&(map_poses[frame_id][0]), R_refined.data());
        t_refined = Vec3(map_poses[frame_id][3], map_poses[frame_id][4], map_poses[frame_id][5]);
        s_refined = map_poses[frame_id][6];

        // Update the pose
        frame->setPose_cr_Rts(R_refined,t_refined,s_refined);
      }
    }
    if (frame_i && !frame_i->isActive())
    {
      // Get frame and update values
      size_t frame_id = frame_i->getFrameId();
      ceres::AngleAxisToRotationMatrix(&(map_poses[frame_id][0]), R_refined.data());
      t_refined = Vec3(map_poses[frame_id][3], map_poses[frame_id][4], map_poses[frame_id][5]);
      s_refined = map_poses[frame_id][6];

      // Update the pose
      frame_i->setPose_cr_Rts(R_refined,t_refined,s_refined);
    }
    return true;
  }
}


}
}
