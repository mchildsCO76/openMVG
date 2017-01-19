
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <openMVG/vsslam/mapping/Cartographer.hpp>

using namespace openMVG;

namespace openMVG  {
namespace VSSLAM  {

// Constructor
Cartographer::Cartographer()
{

}

void Cartographer::addCameraIntrinsicData(const size_t & cam_id, IntrinsicBase * cam_intrinsic)
{
  slam_data.cam_intrinsics[cam_id] = cam_intrinsic;
}

void Cartographer::addObservationToLandmark(const size_t & landmark_id, const std::shared_ptr<Frame> & frame, const size_t featId)
{
  slam_data.structure[landmark_id].obs_[frame.get()] = featId;
}

void Cartographer::addKeyFrameToMap(const std::shared_ptr<Frame> & frame)
{
  slam_data.keyframes[frame->frameId_] = frame->share_ptr();
  // Verbose
  std::cout<<"MAP: Keyframe "<<frame->frameId_<<" inserted!\n";
}
void Cartographer::addLandmarkToMap(MapLandmark & landmark)
{
  // Get new Id for a landmark

  landmark.id_ = slam_data.next_free_landmark_id++;
  slam_data.structure[landmark.id_] = landmark;
  // Verbose
  std::cout<<"MAP: Landmark "<<slam_data.structure.at(slam_data.structure.size()-1).id_<<" inserted!\n";
}

bool Cartographer::initializeMap
(
  std::shared_ptr<Frame> & frame_1,
  std::shared_ptr<Frame> & frame_2,
  Hash_Map<size_t,size_t> & matches_2_1_idx,
  std::vector<size_t> & inliers
)
{
  // -------------------
  // -- Do tiny BA to refine initial poses
  // -------------------


  // Initial VSSLAM_Data system
  VSSLAM_Data init_map;

  // Add initial keyframes
  init_map.keyframes[frame_1->frameId_] = frame_1->share_ptr();
  init_map.keyframes[frame_2->frameId_] = frame_2->share_ptr();

  init_map.cam_intrinsics[frame_1->camId_] = frame_1->cam_intrinsic_;
  init_map.cam_intrinsics[frame_2->camId_] = frame_2->cam_intrinsic_;
  // Add points to the system
  size_t pt_i = 0;
  for (size_t k : inliers)
  {
    Hash_Map<size_t,size_t>::iterator m_iter = matches_2_1_idx.begin();
    std::advance(m_iter,k);

    // 3D point
    MapLandmark l;
    l.id_ = init_map.next_free_landmark_id++;
    // Position of detected features
    const Vec2
      & x1_ = frame_1->getFeaturePosition(m_iter->second),
      & x2_ = frame_2->getFeaturePosition(m_iter->first);

    // Triangulate results
    TriangulateDLT(frame_1->P_, x1_, frame_2->P_, x2_, &l.pt_);

    // Add observations
    l.obs_[frame_1.get()] = m_iter->second;
    l.obs_[frame_2.get()] = m_iter->first;


    init_map.structure[l.id_] = l;
  }

  std::cout<<"Initial map: Keyframes: "<<init_map.keyframes.size()<<" Landmarks: "<<init_map.structure.size()<<"\n";

  // BA - refine only Structure and Rotations & translations (keep intrinsic constant)
  VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
  options.linear_solver_type_ = ceres::DENSE_SCHUR;
  VSSLAM_Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
  if (!bundle_adjustment_obj.Adjust(init_map,
      sfm::Optimize_Options
      (
        cameras::Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
        sfm::Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
        sfm::Structure_Parameter_Type::ADJUST_ALL) // Adjust structure
      )
    )
  {
    return false;
  }


  return false;

}

void Cartographer::clearMap()
{
  // Clear everything
  slam_data.structure.clear();
  slam_data.keyframes.clear();

}





}
}
