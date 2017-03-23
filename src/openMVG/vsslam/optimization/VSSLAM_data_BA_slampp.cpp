
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <type_traits> //std::conditional

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

#include <openMVG/vsslam/optimization/VSSLAM_data_BA_slampp.hpp>

#include <openMVG/vsslam/optimization/slampp_optimizers/SlamPP_Optimizer_Sim3_gXYZ_gXYZ.hpp>


using namespace openMVG;

namespace openMVG {
namespace VSSLAM {



VSSLAM_Bundle_Adjustment_SlamPP::VSSLAM_Bundle_Adjustment_SlamPP
(
  BA_options_SlamPP options
)
:options_(options)
{
  if (options_.map_landmark_type_ == MAP_POINT_TYPE::GLOBAL_EUCLIDEAN && options_.map_camera_type_ == MAP_CAMERA_TYPE::GLOBAL)
  {
    problem_ = std::unique_ptr<SlamPP_Optimizer>(new SlamPP_Optimizer_Sim3_gXYZ_gXYZ(options_.undefined_cam_id, options_.b_verbose_, options_.b_use_schur_, options_.b_do_marginals_, options_.b_do_icra_style_marginals_));
  }
}

VSSLAM_Bundle_Adjustment_SlamPP::BA_options_SlamPP &
VSSLAM_Bundle_Adjustment_SlamPP::getOptions()
{
  return options_;
}

bool VSSLAM_Bundle_Adjustment_SlamPP::OptimizePose
(
  Frame * frame_i,
  Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_i_idx
)
{
  std::cout<<"Optimize pose: "<<frame_i->getFrameId()<<"\n";
/*
  // if vector of frames is empty and no frame_i is given -> no pose to optimize
  if (!frame_i || matches_3D_pts_frame_i_idx.empty())
  {
    return true;
  }

  // Create a local problem
  size_t next_vertex_idx_slampp = 0;
  size_t next_const_vertex_idx_slampp = std::numeric_limits<size_t>::max(); // Constant vertices have decreasing ids starting with max

  BA_options_SlamPP ba_options;
  std::unique_ptr<SlamPP_Optimizer> optimizer = std::unique_ptr<SlamPP_Optimizer>(new SlamPP_Optimizer_Sim3_gXYZ_gXYZ(ba_options.undefined_cam_id, ba_options.b_verbose_, ba_options.b_use_schur_, ba_options.b_do_marginals_, ba_options.b_do_icra_style_marginals_));

  // Set initial settings
  optimizer->Set_AllBatch(ba_options.b_all_batch);
  optimizer->Set_UpdateThreshold(ba_options.f_update_thresh);
  optimizer->Set_TrustRadius(ba_options.f_trust_radius);
  optimizer->Set_TrustRadius_Persistence(ba_options.b_trust_radius_persistent);


  std::cout<<"Optimize a pose settings\n";
  Hash_Map<IndexT, IndexT> map_frame_id_slampp_id;
  Hash_Map<IndexT, IndexT> map_landmark_id_slampp_id;

  // Add frame of interest (current frame)
  // Matches 3D-2D are relating to this frame

  const IndexT & frame_i_id = frame_i->getFrameId();
  const IndexT & cam_i_idx = frame_i->getCamId();
  IntrinsicBase * & cam_i_intrinsic = frame_i->getCameraIntrinsics();

  // Pose parameters
  Eigen::Matrix<double, 12, 1> frame_i_state;
  frame_i->getPoseInverse_StateVector(frame_i_state);

  std::cout<<"B Frame: "<<frame_i->getFrameId()<<"\n";
  std::cout<<"B T: "<<frame_i->T_cr_<<"\n";
  std::cout<<"B State: "<<frame_i_state<<"\n";

  // Add current camera to local BA
  size_t pose_cam_i_id_slampp = next_vertex_idx_slampp;
  double * frame_i_state_ptr = optimizer->Add_CamVertex(pose_cam_i_id_slampp,frame_i_state);
  next_vertex_idx_slampp++;

  std::cout<<"Optimize a pose cam\n";

  // Add all 2D-3D relations of the free pose
  for (auto & match : matches_3D_pts_frame_i_idx)
  {
    // add landmarks as global point : undefined cam id as landmark has no owner no owner
    optimizer->Add_XYZVertexFixed(next_const_vertex_idx_slampp,ba_options.undefined_cam_id, match.first->X_);

    // add reprojection edge
    optimizer->Add_P2CSim3GEdge(next_const_vertex_idx_slampp,pose_cam_i_id_slampp,frame_i->getFeaturePosition(match.second), frame_i->getFeatureInformationMatrix(match.second));

    // Decrease the const id
    next_const_vertex_idx_slampp--;
  }

  // Optimize the solution
  optimizer->Optimize(ba_options.n_max_inc_iters,ba_options.f_inc_nlsolve_thresh, ba_options.f_inc_nlsolve_thresh);

  // Update back the result of the optimization
  Eigen::Map<const Eigen::VectorXd> frame_i_state_updated(frame_i_state_ptr, 7);
  std::cout<<"B State After: "<<frame_i_state_updated<<"\n";

  frame_i->setPose_rc_sim3(frame_i_state_updated);
  std::cout<<"B T after: "<<frame_i->T_cr_<<"\n";
*/
  return true;
}

bool VSSLAM_Bundle_Adjustment_SlamPP::OptimizeLocal
(
  Hash_Map<Frame*, size_t> & tmp_frames,
  Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> > & tmp_structure,
  Frame * frame_i,
  std::vector<std::unique_ptr<MapLandmark> > & vec_triangulated_pts
)
{
 /* std::cout<<"Optimize Local\n";
  BA_options_SlamPP ba_options;
  // if vector of frames is empty and no frame_i is given -> no pose to optimize
  if (!frame_i)
  {
    return true;
  }
  // No structure to optimize
  if (tmp_structure.empty() && vec_triangulated_pts.empty())
  {
    return true;
  }

  std::unique_ptr<SlamPP_Optimizer> optimizer = std::unique_ptr<SlamPP_Optimizer>(new SlamPP_Optimizer_Sim3_gXYZ_gXYZ(ba_options.b_verbose_, ba_options.b_use_schur_, ba_options.b_do_marginals_, ba_options.b_do_icra_style_marginals_));

  optimizer->Set_AllBatch(ba_options.b_all_batch);
  optimizer->Set_UpdateThreshold(ba_options.f_update_thresh);
  optimizer->Set_TrustRadius(ba_options.f_trust_radius);
  optimizer->Set_TrustRadius_Persistence(ba_options.b_trust_radius_persistent);


  size_t slampp_id = 0;
  size_t const_slampp_id = std::numeric_limits<size_t>::max()-1;

  Hash_Map<Frame*, double*> map_frame_slampp;
  Hash_Map<MapLandmark*, double*> map_landmark_slampp;

  Hash_Map<Frame*, IndexT> map_frame_slampp_id;
  Hash_Map<MapLandmark*, IndexT> map_landmark_slampp_id;

  Eigen::Matrix<double, 12, 1> cam_state;
  for (auto & it_frame : tmp_frames)
  {
    Frame * frame = it_frame.first;
    const IndexT & frameId = frame->getFrameId();
    const IndexT & cam_idx = frame->getCamId();

    // Pose parameters
    frame->getStateVector_rc_sim3(cam_state);


    std::cout<<"C Frame: "<<frame_i->getFrameId()<<"\n";
    std::cout<<"C State: "<<cam_state<<"\n";


    // Add camera to BA (if active it is fixed)
    if (frame->isActive())
    {
      // Negative ids mean
      double * state_ptr = optimizer->Add_CamVertex(const_slampp_id,cam_state);
      map_frame_slampp_id[frame] = const_slampp_id;
      const_slampp_id--;
    }
    else
    {
      double * state_ptr = optimizer->Add_CamVertex(slampp_id,cam_state);
      map_frame_slampp[frame] = state_ptr;
      map_frame_slampp_id[frame] = slampp_id;
      slampp_id++;
    }
  }

  // frame_i
  const IndexT & frame_i_id = frame_i->getFrameId();
  const IndexT & cam_i_idx = frame_i->getCamId();
  IntrinsicBase * & cam_i_intrinsic = frame_i->getCameraIntrinsics();

  // If pose is not already among the others added
  if (map_frame_slampp.find(frame_i) == map_frame_slampp.end())
  {
    // Pose parameters
    frame_i->getStateVector_rc_sim3(cam_state);

    // Add camera to BA (if active it is fixed)
    if (frame_i->isActive())
    {
      // Negative ids mean
      double * state_ptr = optimizer->Add_CamVertex(const_slampp_id,cam_state);
      map_frame_slampp_id[frame_i] = const_slampp_id;
      const_slampp_id--;
    }
    else
    {
      double * state_ptr = optimizer->Add_CamVertex(slampp_id,cam_state);
      map_frame_slampp[frame_i] = state_ptr;
      map_frame_slampp_id[frame_i] = slampp_id;
      slampp_id++;
    }
  }



  // Impose restrictions of local map
  if(!tmp_structure.empty())
  {
    // Add all 2D-3D relations of the free pose
    for (auto & t_landmark : tmp_structure)
    {
      // Build the residual block corresponding to the track observation:
      MapLandmark * map_point = t_landmark.first;

      // add landmarks as global point : -1 as no owner
      double * l_ptr = optimizer->Add_XYZVertex(slampp_id,-1, map_point->X_);
      map_landmark_slampp[map_point] = l_ptr;
      map_landmark_slampp_id[map_point] = slampp_id;
      // Add all observations to
      for (auto & obs : map_point->obs_)
      {
        Frame * & frame_obs = obs.second.frame_ptr;
        optimizer->Add_P2CSim3GEdge(slampp_id,map_frame_slampp_id[frame_obs],frame_obs->getFeaturePosition(obs.second.feat_id), frame_obs->getFeatureInformationMatrix(obs.second.feat_id));
      }
      slampp_id++;
    }
  }

  const size_t frame_i_slampp_id = map_frame_slampp_id[frame_i];
  // Impose 3D-2D matches that are already associated with frame_i (valid points are static)
  for (size_t mp_i = 0; mp_i < frame_i->map_points_.size(); ++mp_i)
  {
    MapLandmark* & map_point = frame_i->map_points_[mp_i];
    if (!map_point || !map_point->isActive())
      continue;

    optimizer->Add_XYZVertexFixed(const_slampp_id,-1, map_point->X_);

    // Build the residual block corresponding to the track observation:
    optimizer->Add_P2CSim3GEdge(const_slampp_id,frame_i_slampp_id,frame_i->getFeaturePosition(mp_i), frame_i->getFeatureInformationMatrix(mp_i));

    const_slampp_id--;
  }

  // Add newly triangulated points (3d points are free)
  if(!vec_triangulated_pts.empty())
  {
    for (std::unique_ptr<MapLandmark> & point_it : vec_triangulated_pts)
    {
      // add landmarks as global point : -1 as no owner
      double * l_ptr = optimizer->Add_XYZVertex(slampp_id,-1, point_it->X_);
      map_landmark_slampp[point_it.get()] = l_ptr;
      map_landmark_slampp_id[point_it.get()] = slampp_id;
      const size_t l_slampp_id = slampp_id;
      slampp_id++;

      LandmarkObservations & obs = point_it->obs_;
      for(auto & m_o : obs)
      {
        Frame * frame_i = m_o.second.frame_ptr;
        const IndexT & frame_i_id = frame_i->getFrameId();

        // Check if frame is added to the system
        // If pose is not already among the others added
        if (map_frame_slampp.find(frame_i) == map_frame_slampp.end())
        {
          // Pose parameters
          frame_i->getStateVector_rc_sim3(cam_state);

          // Add camera to BA (if active it is fixed)
          if (frame_i->isActive())
          {
            // Negative ids mean
            double * state_ptr = optimizer->Add_CamVertex(const_slampp_id,cam_state);
            map_frame_slampp_id[frame_i] = const_slampp_id;
            const_slampp_id--;
          }
          else
          {
            double * state_ptr = optimizer->Add_CamVertex(slampp_id,cam_state);
            map_frame_slampp[frame_i] = state_ptr;
            map_frame_slampp_id[frame_i] = slampp_id;
            slampp_id++;
          }
        }

        // Build the residual block corresponding to the track observation:
        optimizer->Add_P2CSim3GEdge(l_slampp_id,map_frame_slampp_id[frame_i],frame_i->getFeaturePosition(m_o.second.feat_id), frame_i->getFeatureInformationMatrix(m_o.second.feat_id));
      }
    }
  }


  // Optimize the solution
  //optimizer->Optimize(ba_options.n_max_inc_iters,ba_options.f_inc_nlsolve_thresh, ba_options.f_inc_nlsolve_thresh);


  // Update camera poses with refined data

  for (auto & it_frame : map_frame_slampp)
  {
    Frame * frame = it_frame.first;

    Eigen::Map<const Eigen::VectorXd> frame_state_after(it_frame.second, 7);

    std::cout<<"Frame: "<<frame->getFrameId()<<"\n";
    std::cout<<"State: "<<frame_state_after<<"\n";
  }

  for (auto & it_landmark : map_landmark_slampp)
  {
    MapLandmark * landmark = it_landmark.first;

    Eigen::Map<const Eigen::VectorXd> landmark_state_after(it_landmark.second, 3);

    std::cout<<"Landmark: "<<landmark->X_<<"\n";
    std::cout<<"State: "<<landmark_state_after<<"\n";
  }

  /*
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

  // Add all poses in local scope
  for (auto & it_frame : tmp_frames)
  {
    Frame * frame = it_frame.first;
    const IndexT & frameId = frame->getFrameId();
    const IndexT & cam_idx = frame->getCamId();

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

  // frame_i
  const IndexT & frame_i_id = frame_i->getFrameId();
  const IndexT & cam_i_idx = frame_i->getCamId();
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
      parameter_block = &map_intrinsics[cam_i_idx][0];
      problem.AddParameterBlock(parameter_block, map_intrinsics[cam_i_idx].size());
      problem.SetParameterBlockConstant(parameter_block);
    }
  }


  // Impose restrictions of local map
  if(!tmp_structure.empty())
  {
    // Add all 2D-3D relations of the free pose
    for (auto & t_landmark : tmp_structure)
    {
      //Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> >

      MapLandmark * map_point = t_landmark.first;
      // Build the residual block corresponding to the track observation:

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      for (auto & obs : map_point->obs_)
      {
        Frame * & frame_obs = obs.second.frame_ptr;
        ceres::CostFunction* cost_function = IntrinsicsToCostFunction(frame_obs->getCameraIntrinsics(), frame_obs->getFeaturePosition(obs.second.feat_id));
        if (cost_function)
        {
          problem.AddResidualBlock(cost_function,
            p_LossFunction,
            &map_intrinsics[frame_obs->getCamId()][0],
            &map_poses[frame_obs->getFrameId()][0],
            map_point->X_.data());
        }
      }
    }
  }


  // Impose 3D-2D matches that are already associated with frame_i (valid points are static)
  for (size_t mp_i = 0; mp_i < frame_i->map_points_.size(); ++mp_i)
  {
    MapLandmark* & map_point = frame_i->map_points_[mp_i];
    if (!map_point || !map_point->isActive())
      continue;

    // Build the residual block corresponding to the track observation:

    // Each Residual block takes a point and a camera as input and outputs a 2
    // dimensional residual. Internally, the cost function stores the observed
    // image location and compares the reprojection against the observation.
    ceres::CostFunction * cost_function = IntrinsicsToCostFunction(cam_i_intrinsic, frame_i->getFeaturePosition(mp_i));


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

  // Add newly triangulated points (3d points are free)
  if(!vec_triangulated_pts.empty())
  {
    for (std::unique_ptr<MapLandmark> & point_it : vec_triangulated_pts)
    {
      LandmarkObservations & obs = point_it->obs_;
      for(auto & m_o : obs)
      {
        Frame * frame_i = m_o.second.frame_ptr;
        const IndexT & frame_i_id = frame_i->getFrameId();

        // Check if frame is added to the system
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
            parameter_block = &map_intrinsics[cam_i_idx][0];
            problem.AddParameterBlock(parameter_block, map_intrinsics[cam_i_idx].size());
            problem.SetParameterBlockConstant(parameter_block);
          }
        }

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
*/





  return true;
}



}
}
