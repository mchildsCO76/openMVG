
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

void Cartographer::addObservations(Frame * frame)
{
  const size_t frame_id = frame->frameId_;
  for (size_t feat_id = 0; feat_id < frame->map_points_.size(); ++feat_id)
  {
    MapLandmark* & map_point = frame->map_points_[feat_id];
    if (!map_point)
      continue;

    MapObservation m1(feat_id,frame);
    m1.pt_ptr = &(frame->getFeaturePositionUndistorted(feat_id));
    map_point->obs_[frame_id] = m1;

    // Update descriptor to the latest
    feature_extractor_->getDescriptorRaw(frame->regions_.get(),feat_id,&map_point->bestDesc_);

  }

  /*for (auto measurement : new_measuremetns)
  {
    MapLandmark * m_landmark = measurement.first;
    const size_t feat_id = measurement.second;
    MapObservation m1(feat_id,frame);
    // we input distorted/undistorted coordinates based on the camera calibration
    m1.pt_ptr = &(frame->getFeaturePositionUndistorted(feat_id));
    m_landmark->obs_[frame_id] = m1;

    // Update descriptor to the latest
    feature_extractor_->getDescriptorRaw(frame->regions_.get(),feat_id,&m_landmark->bestDesc_);

    // Mark the point in frame as tracked
    frame->map_points_[feat_id] = m_landmark;
  }*/
}
void Cartographer::addObservationToLandmark(const size_t & landmark_id, const std::shared_ptr<Frame> & frame, const size_t featId)
{
  //slam_data.structure[landmark_id].obs_[frame->frameId_] = MapObservation(featId,frame.get(),&(frame->getFeaturePosition(featId)));
}

void Cartographer::addKeyFrameToMap(const std::shared_ptr<Frame> & frame)
{
  slam_data.keyframes[frame->frameId_] = frame->share_ptr();
  slam_data.cam_intrinsics[frame->cam_->cam_id] = frame->cam_->cam_intrinsic_ptr;
  // Verbose
  //std::cout<<"MAP: Keyframe "<<frame->frameId_<<" inserted!\n";
}
MapLandmark * Cartographer::addLandmarkToMap(MapLandmark & landmark)
{
  // Get new Id for a landmark

  landmark.id_ = slam_data.getNextFreeLandmarkId();
  slam_data.structure[landmark.id_] = landmark;
  // Verbose
  //std::cout<<"MAP: Landmark "<<slam_data.structure.at(slam_data.structure.size()-1).id_<<" inserted!\n";
  return & slam_data.structure[landmark.id_];
}

void Cartographer::initializeMap
(
  std::shared_ptr<Frame> & frame_1,
  std::shared_ptr<Frame> & frame_2,
  std::vector<std::pair<Vec3, std::deque<std::pair<Frame*,size_t> > > > & vec_new_pts_3D_obs
)
{
  const size_t frame_id_1 = frame_1->frameId_;
  const size_t frame_id_2 = frame_2->frameId_;
  // Add initial frames to the map
  addKeyFrameToMap(frame_1);
  addKeyFrameToMap(frame_2);
  // Add each successfully triangulated point
  for ( auto new_pt : vec_new_pts_3D_obs)
  {
    MapLandmark landmark;
    MapLandmark * m_landmark = addLandmarkToMap(landmark);
    m_landmark->X_ =  new_pt.first;

    std::deque<std::pair<Frame*,size_t> > & measurements = new_pt.second;
    for (auto measurement: measurements)
    {
      Frame * frame = measurement.first;
      const size_t frameId = frame->frameId_;
      const size_t & feat_id = measurement.second;

      MapObservation m1(feat_id,frame);
      m1.pt_ptr = &(frame->getFeaturePositionUndistorted(feat_id));
      m_landmark->obs_[frameId] = m1;
      frame->map_points_[feat_id] = m_landmark;
    }

    std::pair<Frame*,size_t> & last_measure = new_pt.second.back();
    
    feature_extractor_->getDescriptorRaw(last_measure.first->regions_.get(),last_measure.second,&m_landmark->bestDesc_);
  }

  // Update covisibility structure between both initial frames
  frame_2->updateFrameVisibilityConnections();

}

void Cartographer::clearMap()
{
  // Clear everything
  slam_data.structure.clear();
  // Clear all pointers in keyframes to points in the map
  for (auto kf: slam_data.keyframes)
  {
    kf.second->clearMapPoints();
  }
  slam_data.keyframes.clear();
  slam_data.next_free_landmark_id = 0;

}


void Cartographer::updateLocalMap(Frame * currentFrame)
{
  local_map_points_.clear();
  // Get neighbour frames of the current frame
  local_map_frames_ = currentFrame->getBestCovisibilityFrames(10);

  // Get all already reconstructed points from all the frames
  for(Frame* neigh_frame : local_map_frames_)
  {
    const size_t frameId = currentFrame->frameId_;
    for(MapLandmark * map_point: neigh_frame->map_points_)
    {
      if (!map_point || map_point->localMapFrameId_ == frameId)
        continue;
      local_map_points_.push_back(map_point);
      map_point->localMapFrameId_ = frameId;
    }
  }
}

std::vector<MapLandmark*> Cartographer::getLocalMapPoints(Frame * currentFrame, std::vector<Frame*> & neighbor_frames)
{
  std::vector<MapLandmark*> local_map_points;

  // Get all already reconstructed points from all the frames
  for(Frame* neigh_frame : neighbor_frames)
  {
    const size_t frameId = currentFrame->frameId_;
    for(MapLandmark * map_point: neigh_frame->map_points_)
    {
      if (!map_point || map_point->localMapFrameId_ == frameId)
        continue;
      local_map_points.push_back(map_point);
      map_point->localMapFrameId_ = frameId;
    }
  }
  return local_map_points;
}




}
}
