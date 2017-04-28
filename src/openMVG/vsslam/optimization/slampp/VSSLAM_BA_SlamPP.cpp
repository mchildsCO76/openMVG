// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <openMVG/vsslam/optimization/slampp/VSSLAM_BA_SlamPP.hpp>


namespace openMVG {
namespace vsslam {



VSSLAM_BA_SlamPP::VSSLAM_BA_SlamPP
(
  BA_options_SlamPP options
)
:options_(options)
{
  //if (options_.b_use_loss_function_)

  problem_ = std::unique_ptr<SlamPP_Optimizer>(new SlamPP_Optimizer_Sim3_gXYZ_gXYZ(options_.undefined_cam_id, options_.b_verbose_, options_.b_use_schur_, options_.b_do_marginals_, options_.b_do_icra_style_marginals_));
  // Set initial settings
  problem_->Set_AllBatch(options_.b_all_batch);
  problem_->Set_UpdateThreshold(options_.f_update_thresh);
  problem_->Set_TrustRadius(options_.f_trust_radius);
  problem_->Set_TrustRadius_Persistence(options_.b_trust_radius_persistent);

}

VSSLAM_BA_SlamPP::BA_options_SlamPP & VSSLAM_BA_SlamPP::getOptions()
{
  return options_;
}

// Local optimization
bool VSSLAM_BA_SlamPP::OptimizeLocalSystem
(
  Frame * frame_i,
  NewMapLandmarks & vec_new_landmarks,
  bool b_use_loss_function,
  BA_options_SlamPP & ba_options
)
{

  return false;
}

// Pose optimization
bool VSSLAM_BA_SlamPP::OptimizePose
(
  Frame * frame,
  Hash_Map<MapLandmark *,IndexT> & matches_map_cur_idx,
  bool b_use_loss_function,
  BA_options_SlamPP & ba_options
)
{
  return false;
}


bool VSSLAM_BA_SlamPP::addFrameToGlobalSystem(Frame * frame, bool b_frame_fixed)
{

  // Add frame of interest
  const IndexT & frame_id = frame->getFrameId();
  const IndexT & cam_id = frame->getCamId();

  if (map_poses_.find(frame_id)!=map_poses_.end())
  {
    std::cout<<"Cartographer: [Slam++ GlobalBA] Frame "<<frame->getFrameId()<< "already in the system!";
    return false;
  }

  Eigen::Matrix<double, 12, 1> vec_state_frame;
  // Get pose in the WORLD reference frame
  frame->getPose_StateVector(vec_state_frame,nullptr);

  size_t frame_slampp_id = getNextVertexId();
  //std::cout<<"Slam++: Add frame id: "<<frame_slampp_id<<"\n";
  double * ptr_state_frame;
  if (b_frame_fixed)
  {
    // Set the observing camera as fixed
    ptr_state_frame = problem_->Add_CamVertexFixed(frame_slampp_id,vec_state_frame);
  }
  else
  {
    ptr_state_frame = problem_->Add_CamVertex(frame_slampp_id,vec_state_frame);
  }
  // Add camera to map
      map_poses_[frame_id] = std::make_pair(frame_slampp_id,ptr_state_frame);

  std::cout<<"Cartographer: [Slam++ GlobalBA] Add frame: "<<frame->getFrameId()<< " Fixed: "<<b_frame_fixed<<" to global map!\n";
  return true;

}

bool VSSLAM_BA_SlamPP::addLandmarkToGlobalSysyem(MapLandmark * map_landmark)
{
  // Create vertex for Landmark
  const size_t landmark_slampp_id = getNextVertexId();

  //std::cout<<"Slam++: Add landmark id: "<<landmark_slampp_id<<"\n";
  // Add landmarks as global point : m_undefined_camera_id as no owner
  double * landmark_ptr = problem_->Add_XYZVertex(landmark_slampp_id,options_.undefined_cam_id, map_landmark->X_);

  // Add landmark to map
  map_landmarks_[map_landmark->id_] = std::make_pair(landmark_slampp_id,landmark_ptr);


  LandmarkObservations & vec_obs = map_landmark->getObservations();
  for(auto & m_o : vec_obs)
  {
    Frame * frame = m_o.second.frame_ptr;
    const IndexT & feat_id_frame = m_o.second.feat_id;

    const IndexT & frame_id = frame->getFrameId();
    const IndexT & cam_id_frame = frame->getCamId();
    IntrinsicBase * & cam_intrinsic = frame->getCameraIntrinsics();

    // Add frame to the problem if its not added yet
    if (map_poses_.find(frame_id) == map_poses_.end())
    {
      std::cout<<"Cartographer: [SlamPP GlobalBA] Adding landmark: Frame "<<frame_id<<" is not yet in the system!! Skipping landmark!";
      return false;
    }

    // Add measurement edge
    Eigen::Matrix2d inf_mat = frame->getFeatureSqrtInfMatrix(feat_id_frame).cwiseProduct(frame->getFeatureSqrtInfMatrix(feat_id_frame));
    const size_t frame_slampp_id = map_poses_.find(frame_id)->second.first;
    problem_->Add_P2CSim3GEdge(landmark_slampp_id,frame_slampp_id,frame->getFeaturePosition(feat_id_frame), inf_mat);

  }
  return true;
}

bool VSSLAM_BA_SlamPP::addObservationToGlobalSystem(MapLandmark * map_landmark, MapObservation * map_observation)
{

  Frame * frame = map_observation->frame_ptr;
  const IndexT & feat_id_frame = map_observation->feat_id;

  const IndexT & frame_id = frame->getFrameId();
  const IndexT & landmark_id = map_landmark->id_;
  const IndexT & cam_id_frame = frame->getCamId();
  IntrinsicBase * & cam_intrinsic = frame->getCameraIntrinsics();

  // Add frame to the problem if its not added yet
  if (map_poses_.find(frame_id) == map_poses_.end())
  {
    std::cout<<"Cartographer: [SlamPP GlobalBA] Adding landmark: Frame "<<frame_id<<" is not yet in the system!! Skipping landmark!";
    return false;
  }

  // Check if landmark exists in the system
  if (map_landmarks_.find(landmark_id) == map_landmarks_.end() || !map_landmark->isActive())
  {
    std::cout<<"Cartographer: [SlamPP GlobalBA] Adding landmark: "<<landmark_id<<" is not yet in the system!! Skipping landmark!";
    return false;
  }

  if (map_observations_.find(frame_id)!=map_observations_.end())
  {
    bool b_exists = false;
    for (auto mo : map_observations_[frame_id])
    {
      if (mo == landmark_id)
      {
        b_exists = true;
        break;
      }
    }

    if (b_exists)
    {
      std::cout<<"OBSERVATION EXISTS!!!!\n";
      return false;
    }
    else
    {
      map_observations_[frame_id].push_back(landmark_id);
    }
  }
  else
  {
    map_observations_[frame_id].push_back(landmark_id);
  }

  // Add measurement edge
  const size_t frame_slampp_id = map_poses_.find(frame_id)->second.first;
  const size_t landmark_slampp_id = map_landmarks_.find(landmark_id)->second.first;


  //std::cout<<"Slam++: Add observation: landmark id: "<<landmark_slampp_id<<" frame: "<<frame_slampp_id<<"\n";
  // Add measurement edge
  Eigen::Matrix2d inf_mat = frame->getFeatureSqrtInfMatrix(feat_id_frame).cwiseProduct(frame->getFeatureSqrtInfMatrix(feat_id_frame));
  problem_->Add_P2CSim3GEdge(landmark_slampp_id,frame_slampp_id,frame->getFeaturePosition(feat_id_frame), inf_mat);

  return true;
}

bool VSSLAM_BA_SlamPP::optimizeGlobal(VSSLAM_Map & map_global)
{
  for (auto & frame_it : map_poses_)
  {
    const IndexT & frame_id = frame_it.first;
    Frame * frame = map_global.map_frame_[frame_id].get();
    std::cout<<"Frame before: "<<frame->getFrameId()<<" T: "<<frame->getTransformationMatrix()<<"\n";

    Eigen::Matrix<double, 12, 1> vec_state_frame;
    // Get pose in the WORLD reference frame
    frame->getPose_StateVector(vec_state_frame,nullptr);
    std::cout<<"sim3: "<<vec_state_frame<<"\n";
    std::cout<<"Scale: "<<frame->getPoseScale()<<"\n";
  }

  // Optimize the solution
  problem_->Optimize(options_.n_max_inc_iters,options_.f_inc_nlsolve_thresh, options_.f_inc_nlsolve_thresh);

  std::cout<<"UPDATE DATA\n";

  // Update frames with refined data
  for (auto & frame_it : map_poses_)
  {
    const IndexT & frame_id = frame_it.first;
    Frame * frame = map_global.map_frame_[frame_id].get();

    Eigen::Map<Eigen::VectorXd> frame_state_after = problem_->r_Vertex_State(frame_it.second.first);;
    Eigen::Matrix<double, 7, 1>  vec_state = frame_state_after.template head<7>();
    frame->setPose_sim3(vec_state,frame->getReferenceFrame());
    std::cout<<"Frame: "<<frame->getFrameId()<<" T: "<<frame->getTransformationMatrix()<<"\n";
    std::cout<<"sim3: "<<vec_state<<"\n";
    std::cout<<"Scale: "<<frame->getPoseScale()<<"\n";
  }

  // Upate landmarks with refined data

  for (auto & landmark_it : map_landmarks_)
  {
    const IndexT & landmark_id = landmark_it.first;
    const size_t & landmark_slampp_id = landmark_it.second.first;
    // Get landmark pointer from the structure
    MapLandmark * landmark = map_global.getLandmark(landmark_id).get();


    // Recover value from slampp
    Eigen::Map<Eigen::VectorXd> landmark_state_after = problem_->r_Vertex_State(landmark_slampp_id);
    landmark->X_ = landmark_state_after;
    //std::cout<<"Landmark id:"<<landmark->id_<<" X: "<<landmark->X_<<"\n";
  }

  std::cout<<"Cartographer: [Slam++ GlobalBA] Optimized OK!\n";


  return true;
}

}
}
