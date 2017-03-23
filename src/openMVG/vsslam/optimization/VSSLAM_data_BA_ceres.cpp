
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



#include <openMVG/cameras/Camera_Common.hpp>

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
#include <openMVG/vsslam/optimization/sim3.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA_ceres.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA_ceres_camera_functor.hpp>
//#include <openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp>
using namespace openMVG;

namespace openMVG {
namespace VSSLAM {

/// Create the appropriate cost functor according the provided input camera intrinsic model.
/// The residual can be weighetd if desired (default 0.0 means no weight).
ceres::CostFunction * IntrinsicsToCostFunction
(
  IntrinsicBase * intrinsic,
  const Vec2 & observation,
  const Eigen::Matrix<double, 2, 2> & inf_matrix,
  const double weight
)
{
  switch(intrinsic->getType())
  {
    case PINHOLE_CAMERA:
      //return ResidualErrorFunctor_Pinhole_Intrinsic_Rts::Create(observation, weight);
      return Chi2ErrorFunctor_Pinhole_Intrinsic_Rts::Create(observation, inf_matrix, weight);
     break;
    default:
      return nullptr;
  }
}

VSSLAM_Bundle_Adjustment_Ceres::VSSLAM_Bundle_Adjustment_Ceres
(
  BA_options_Ceres options
)
:options_(options)
{
}

VSSLAM_Bundle_Adjustment_Ceres::BA_options_Ceres &
VSSLAM_Bundle_Adjustment_Ceres::getOptions()
{
  return options_;
}


bool VSSLAM_Bundle_Adjustment_Ceres::OptimizePose
(
  Frame * frame_i,
  Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_i_idx
)
{
  std::cout<<"Optimize pose: "<<frame_i->getFrameId()<<"\n";

  // if vector of frames is empty and no frame_i is given -> no pose to optimize
  if (!frame_i || matches_3D_pts_frame_i_idx.empty())
  {
    return true;
  }

  // Create the local problem
  BA_options_Ceres ba_options;
  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));

  // Pose parameters
  Mat3 R;
  double angleAxis[3];
  Vec3 t;
  double s;

  // Add frame of interest (current frame)
  // Matches 3D-2D are relating to this frame
  const IndexT & frame_i_id = frame_i->getFrameId();
  const IndexT & cam_i_idx = frame_i->getCamId();
  IntrinsicBase * & cam_i_intrinsic = frame_i->getCameraIntrinsics();

  // Get pose in the WORLD reference frame
  frame_i->getPose_Rts(R,t,s,nullptr);

  std::cout<<"Pose before: \n R: "<<R<<"\n t: "<<t<<"\n s: "<<s<<"\n";

  ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
  // Save data about the poses
  map_poses[frame_i_id] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2),s};
  // Add parameter to the block
  double * parameter_block = &map_poses[frame_i_id][0];
  problem.AddParameterBlock(parameter_block, 7);
  // Add intrinsic parameters
  map_intrinsics[cam_i_idx] = frame_i->getCameraIntrinsics()->getParams();
  parameter_block = &map_intrinsics[cam_i_idx][0];
  problem.AddParameterBlock(parameter_block, map_intrinsics[cam_i_idx].size());
  // If camera is calibrated we fix the parameters
  if (frame_i->getCamCalibrated())
    problem.SetParameterBlockConstant(parameter_block);

  // Add all 3D-2D matchings -> they are fixed
  for (auto & match : matches_3D_pts_frame_i_idx)
  {
    // Each Residual block takes a point and a camera as input and outputs a 2
    // dimensional residual. Internally, the cost function stores the observed
    // image location and compares the reprojection against the observation.
    ceres::CostFunction* cost_function = IntrinsicsToCostFunction(cam_i_intrinsic,
                                                                  frame_i->getFeaturePosition(match.second),
                                                                  frame_i->getFeatureInformationMatrix(match.second));

    if (cost_function)
    {
      problem.AddResidualBlock(cost_function,
        p_LossFunction,
        &map_intrinsics[cam_i_idx][0],
        &map_poses[frame_i_id][0],
        match.first->X_.data());
      // Set structure fixed if in the global map
      // For now we set everything fixed -> pose of the current frame is just approximated
      // We do this to make the BA really small
      //if (match.first->isActive())
        problem.SetParameterBlockConstant(match.first->X_.data());
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type = ba_options.preconditioner_type_;
  ceres_config_options.linear_solver_type = ba_options.linear_solver_type_;
  ceres_config_options.sparse_linear_algebra_library_type = ba_options.sparse_linear_algebra_library_type_;
  ceres_config_options.minimizer_progress_to_stdout = ba_options.b_verbose_;
  ceres_config_options.logging_type = ceres::SILENT;
  ceres_config_options.num_threads = ba_options.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ba_options.nb_threads_;
  ceres_config_options.parameter_tolerance = ba_options.parameter_tolerance_;

  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ba_options.b_ceres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ba_options.b_verbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ba_options.b_verbose_)
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

    // Update camera poses with refined data
    Mat3 R_refined;
    // Get frame and update values
    IndexT frame_id = frame_i->getFrameId();
    ceres::AngleAxisToRotationMatrix(&(map_poses[frame_id][0]), R_refined.data());
    Vec3 t_refined = Vec3(map_poses[frame_id][3], map_poses[frame_id][4], map_poses[frame_id][5]);
    double s_refined = map_poses[frame_id][6];

    std::cout<<"Pose after: \n R: "<<R_refined<<"\n t: "<<t_refined<<"\n s: "<<s_refined<<"\n";
    // Update the pose
    frame_i->setPose_Rts(R_refined,t_refined,s_refined, nullptr);


    // Get pose in the WORLD reference frame
    frame_i->getPose_Rts(R,t,s,nullptr);
    std::cout<<"Pose after B: \n R: "<<R<<"\n t: "<<t<<"\n s: "<<s<<"\n";
    return true;
  }
}

// OptimizeLocal optimizes a system consisting of:
//  - Pose of current frame
//  - global landmarks (fixed) seen in the frame
//  - local landmarks (free) seen in the frame
//  - any camera pose (fixed) used by local landmarks

bool VSSLAM_Bundle_Adjustment_Ceres::OptimizeLocal
(
  Hash_Map<Frame*, size_t> & tmp_frames,
  Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> > & tmp_structure,
  Frame * frame_i,
  std::vector<std::unique_ptr<MapLandmark> > & vec_triangulated_pts
)
{
  // if vector of frames is empty and no frame_i is given -> no pose to optimize
  if (!frame_i)
  {
    return true;
  }

  std::cout<<"Optimize local map\n";

  // Create the local problem
  BA_options_Ceres ba_options;
  ceres::Problem problem;

  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(4.0));

  // Pose parameters
  Mat3 R;
  double angleAxis[3];
  Vec3 t;
  double s;

  // Add frame of interest (current frame)
  const IndexT & frame_i_id = frame_i->getFrameId();
  const IndexT & cam_i_idx = frame_i->getCamId();
  IntrinsicBase * & cam_i_intrinsic = frame_i->getCameraIntrinsics();

  // Add frame_i to the problem

  // Get pose in the WORLD reference frame
  frame_i->getPose_Rts(R,t,s,nullptr);
  // Convert rotation matrix to angleaxis representation
  ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
  // Save data about the poses
  map_poses[frame_i_id] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2),s};
  // Add parameter to the block
  std::cout<<"Add frame_i POSE\n";
  double * parameter_block = &map_poses[frame_i_id][0];
  problem.AddParameterBlock(parameter_block, 7);
  std::cout<<"Add frame_i POSE A\n";
  // Add intrinsic parameters
  std::cout<<"Add frame_i INTRINSICS\n";
  map_intrinsics[cam_i_idx] = frame_i->getCameraIntrinsics()->getParams();
  parameter_block = &map_intrinsics[cam_i_idx][0];
  problem.AddParameterBlock(parameter_block, map_intrinsics[cam_i_idx].size());
  std::cout<<"Add frame_i INTRINSICS B\n";
  // If camera is calibrated we fix the parameters
  if (frame_i->getCamCalibrated())
    problem.SetParameterBlockConstant(parameter_block);

  std::cout<<"Add frame_i INTRINSICSC\n";


  // Impose 3D-2D matches that are already associated with frame_i
  // Global landmarks are fixed
  // Local landmarks are free -> we add all cameras involved as fixed
  for (size_t mp_i = 0; mp_i < frame_i->map_points_.size(); ++mp_i)
  {
    MapLandmark* & map_point = frame_i->map_points_[mp_i];
    if (!map_point)
      continue;

    // Add measurement in frame_i to landmark
    // Create the cost function for the measurement
    ceres::CostFunction* cost_function = IntrinsicsToCostFunction(cam_i_intrinsic,
        frame_i->getFeaturePosition(mp_i),
        frame_i->getFeatureInformationMatrix(mp_i));

    // Add cost term
    if (cost_function)
    {
      problem.AddResidualBlock(cost_function,
        p_LossFunction,
        &map_intrinsics[cam_i_idx][0],
        &map_poses[frame_i_id][0],
        map_point->X_.data());
      // We do not fix the point as its only in local structure
    }

    if (map_point->isActive())
    {
      // If landmark is global we fix it
      problem.SetParameterBlockConstant(map_point->X_.data());
    }
    else
    {
      // If local map point we add all other cameras that see the point
      for (auto & m_o : map_point->obs_)
      {
        Frame * & frame_obs = m_o.second.frame_ptr;

        // Frame of interest (current frame)
        const IndexT & frame_obs_id = frame_obs->getFrameId();
        const IndexT & cam_obs_idx = frame_obs->getCamId();
        IntrinsicBase * & cam_obs_intrinsic = frame_obs->getCameraIntrinsics();

        // Add frame to the problem if its not added yet
        if (map_poses.find(frame_obs_id) == map_poses.end())
        {
          // Get pose in the WORLD reference frame
          frame_obs->getPose_Rts(R,t,s,nullptr);
          // Convert rotation matrix to angleaxis representation
          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
          // Save data about the poses
          map_poses[frame_obs_id] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2),s};
          // Add parameter to the block
          double * parameter_block = &map_poses[frame_obs_id][0];
          problem.AddParameterBlock(parameter_block, 7);
          // Set the observing camera as fixed
          problem.SetParameterBlockConstant(parameter_block);

          if (map_intrinsics.find(cam_obs_idx) == map_intrinsics.end())
          {
            // Add intrinsic parameters
            map_intrinsics[cam_obs_idx] = frame_obs->getCameraIntrinsics()->getParams();
            parameter_block = &map_intrinsics[cam_obs_idx][0];
            problem.AddParameterBlock(parameter_block, map_intrinsics[cam_obs_idx].size());
            // Fix the calibration parameteres (even if the camera is not calibrated - we are not changing this cam)
            problem.SetParameterBlockConstant(parameter_block);
          }
        }

        // Create the cost function for the measurement
        ceres::CostFunction* cost_function = IntrinsicsToCostFunction(cam_obs_intrinsic,
            frame_obs->getFeaturePosition(m_o.second.feat_id),
            frame_obs->getFeatureInformationMatrix(m_o.second.feat_id));

        // Add cost term
        if (cost_function)
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[cam_obs_idx][0],
            &map_poses[frame_obs_id][0],
            map_point->X_.data());
        }
      }
    }
  }
  std::cout<<"JJJJ\n";

  // Add newly triangulated points (3d points are free)
  // Cameras added are fixed as they are already global
  if(!vec_triangulated_pts.empty())
  {
    for (std::unique_ptr<MapLandmark> & point_it : vec_triangulated_pts)
    {
      std::cout<<"A1\n";
      LandmarkObservations & obs = point_it->obs_;
      for(auto & m_o : obs)
      {

        std::cout<<"A11\n";
        Frame * frame_obs = m_o.second.frame_ptr;

        // Frame of interest (current frame)
        const IndexT & frame_obs_id = frame_obs->getFrameId();
        const IndexT & cam_obs_idx = frame_obs->getCamId();
        IntrinsicBase * & cam_obs_intrinsic = frame_obs->getCameraIntrinsics();

        // Add frame to the problem if its not added yet
        if (map_poses.find(frame_obs_id) == map_poses.end())
        {
          // Get pose in the WORLD reference frame
          frame_obs->getPose_Rts(R,t,s,nullptr);
          // Convert rotation matrix to angleaxis representation
          ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
          // Save data about the poses
          map_poses[frame_obs_id] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2),s};
          // Add parameter to the block
          double * parameter_block = &map_poses[frame_obs_id][0];
          problem.AddParameterBlock(parameter_block, 7);
          // Set the observing camera as fixed
          problem.SetParameterBlockConstant(parameter_block);

          if (map_intrinsics.find(cam_obs_idx) == map_intrinsics.end())
          {
            // Add intrinsic parameters
            map_intrinsics[cam_obs_idx] = frame_obs->getCameraIntrinsics()->getParams();
            parameter_block = &map_intrinsics[cam_obs_idx][0];
            problem.AddParameterBlock(parameter_block, map_intrinsics[cam_obs_idx].size());
            // Fix the calibration parameteres (even if the camera is not calibrated - we are not changing this cam)
            problem.SetParameterBlockConstant(parameter_block);
          }
        }

        std::cout<<"A12\n";
        // Create the cost function for the measurement
        ceres::CostFunction* cost_function = IntrinsicsToCostFunction(cam_obs_intrinsic,
            frame_obs->getFeaturePosition(m_o.second.feat_id),
            frame_obs->getFeatureInformationMatrix(m_o.second.feat_id));

        std::cout<<"A13\n";
        // Add cost term
        if (cost_function)
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[cam_obs_idx][0],
            &map_poses[frame_obs_id][0],
            point_it->X_.data());
        }

        std::cout<<"A14\n";
      }

      std::cout<<"A2\n";
    }
  }
  std::cout<<"KKK\n";


  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type = ba_options.preconditioner_type_;
  ceres_config_options.linear_solver_type = ba_options.linear_solver_type_;
  ceres_config_options.sparse_linear_algebra_library_type = ba_options.sparse_linear_algebra_library_type_;
  ceres_config_options.minimizer_progress_to_stdout = ba_options.b_verbose_;
  ceres_config_options.logging_type = ceres::SILENT;
  ceres_config_options.num_threads = ba_options.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ba_options.nb_threads_;
  ceres_config_options.parameter_tolerance = ba_options.parameter_tolerance_;

  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ba_options.b_ceres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ba_options.b_verbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ba_options.b_verbose_)
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

    // Update camera poses with refined data
    Mat3 R_refined;
    // Get frame and update values
    IndexT frame_id = frame_i->getFrameId();
    ceres::AngleAxisToRotationMatrix(&(map_poses[frame_id][0]), R_refined.data());
    Vec3 t_refined = Vec3(map_poses[frame_id][3], map_poses[frame_id][4], map_poses[frame_id][5]);
    double s_refined = map_poses[frame_id][6];

    std::cout<<"Pose afterFF: \n R: "<<R_refined<<"\n t: "<<t_refined<<"\n s: "<<s_refined<<"\n";
    // Update the pose
    frame_i->setPose_Rts(R_refined,t_refined,s_refined, nullptr);


    // Get pose in the WORLD reference frame
    frame_i->getPose_Rts(R,t,s,nullptr);
    std::cout<<"Pose after FFB: \n R: "<<R<<"\n t: "<<t<<"\n s: "<<s<<"\n";
    return true;
  }

  return true;
}


}
}
