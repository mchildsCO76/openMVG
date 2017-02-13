
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <openMVG/vsslam/mapping/Cartographer.hpp>

using namespace openMVG;

namespace openMVG  {
namespace VSSLAM  {

// Constructor
Cartographer::Cartographer(){}

void Cartographer::AddLandmarksToMap
(
  std::vector<std::pair<Vec3, std::deque<std::pair<Frame*,size_t> > > > & vec_new_pts_3D_obs
)
{
  // Add each successfully triangulated point
  for ( auto new_pt : vec_new_pts_3D_obs)
  {
    MapLandmark * m_landmark = newLandmarkToMap();
    m_landmark->X_ =  new_pt.first;
    //std::cout<<"X: "<<m_landmark->X_<<"\n";

    std::deque<std::pair<Frame*,size_t> > & measurements = new_pt.second;
    for (auto measurement: measurements)
    {
      Frame * frame = measurement.first;
      const size_t & frameId = frame->getFrameId();
      const size_t & feat_id = measurement.second;
      //std::cout<<"pt: "<<frameId<<" :: "<<feat_id<<"\n";
      MapObservation m1(feat_id,frame);
      m1.pt_ptr = &(frame->getFeaturePositionUndistorted(feat_id));
      //std::cout<<"P: "<<*m1.pt_ptr<<"\n";
      m_landmark->obs_[frameId] = m1;
      frame->map_points_[feat_id] = m_landmark;
    }

    std::pair<Frame*,size_t> & last_measure = new_pt.second.back();

    feature_extractor_->getDescriptorRaw(last_measure.first->regions_.get(),last_measure.second,& m_landmark->bestDesc_);
  }
}
void Cartographer::initializeMap
(
  std::shared_ptr<Frame> & frame_1,
  std::shared_ptr<Frame> & frame_2,
  std::vector<std::pair<Vec3, std::deque<std::pair<Frame*,size_t> > > > & vec_new_pts_3D_obs
)
{
  std::cout<<"Cartographer: Initialize map!\n";

  const size_t & frame_id_1 = frame_1->getFrameId();
  const size_t & frame_id_2 = frame_2->getFrameId();

  // Add initial frames to the map
  addKeyframeToMap(frame_1);
  addKeyframeToMap(frame_2);

  // Add each successfully triangulated point
  for ( auto new_pt : vec_new_pts_3D_obs)
  {
    MapLandmark * m_landmark = newLandmarkToMap();
    m_landmark->X_ =  new_pt.first;

    std::deque<std::pair<Frame*,size_t> > & measurements = new_pt.second;
    for (auto measurement: measurements)
    {
      Frame * frame = measurement.first;
      const size_t & frameId = frame->getFrameId();
      const size_t & feat_id = measurement.second;

      MapObservation m1(feat_id,frame);
      m1.pt_ptr = &(frame->getFeaturePositionUndistorted(feat_id));
      m_landmark->obs_[frameId] = m1;
      frame->map_points_[feat_id] = m_landmark;
    }

    std::pair<Frame*,size_t> & last_measure = new_pt.second.back();

    feature_extractor_->getDescriptorRaw(last_measure.first->regions_.get(),last_measure.second,& m_landmark->bestDesc_);
  }

  // Update covisibility structure between both initial frames
  frame_2->computeFrameVisibilityConnections();
}

// Clears all points and frames in the map (keeps camera intrinsics)
void Cartographer::clearMap()
{
  // Clear everything
  structure.clear();
  // Clear all pointers in keyframes to points in the map
  for (auto kf: keyframes)
  {
    kf.second->clearMapPoints();
  }
  keyframes.clear();
  next_free_landmark_id = 0;

}


// ------------------------------
// -- Add data to map
// ------------------------------

// Add camera to map data
void Cartographer::addCameraData(const size_t & cam_id, IntrinsicBase * cam_intrinsic)
{
  cam_intrinsics[cam_id] = cam_intrinsic;
}

void Cartographer::addKeyframeToMap(const std::shared_ptr<Frame> & frame)
{
  std::cout<<"Cartographer: Add keyframe "<<frame->getFrameId()<<" to map!\n";
  keyframes[frame->getFrameId()] = frame->share_ptr();
  cam_intrinsics[frame->getCamId()] = frame->getCameraIntrinsics();
}

MapLandmark * Cartographer::newLandmarkToMap()
{
  // Get new Id for a landmark
  /*
  MapLandmark landmark;
  landmark.id_ = getNextFreeLandmarkId();
  structure[landmark.id_] = landmark;*/
  // Verbose
  //std::cout<<"MAP: Landmark "<<slam_data.structure.at(slam_data.structure.size()-1).id_<<" inserted!\n";
  return & structure[getNextFreeLandmarkId()];
}

void Cartographer::addObservationsToLandmarks(Frame * frame)
{
  std::cout<<"Cartographer: Add observations to landmarks! From frame: "<<frame->getFrameId()<<"!\n";
  const size_t & frame_id = frame->getFrameId();
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
}


// ------------------------------
// -- Local map
// ------------------------------
void Cartographer::getLocalMapPoints(Frame * currentFrame, std::vector<Frame*> & local_frames, std::vector<MapLandmark*> & local_points)
{
  const size_t & frame_cur_id = currentFrame->getFrameId();

  // Get all already reconstructed points from all the frames
  for(Frame* neigh_frame : local_frames)
  {
    for(MapLandmark * map_point: neigh_frame->map_points_)
    {
      // If triangulated and not yet accounted for
      if (!map_point || map_point->localMapFrameId_ == frame_cur_id)
        continue;
      local_points.push_back(map_point);
      map_point->localMapFrameId_ = frame_cur_id;
    }
  }
}

/*
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
*/



}
}
