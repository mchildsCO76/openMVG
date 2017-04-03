
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/vsslam/Frame.hpp>

namespace openMVG  {
namespace VSSLAM  {

class Frame;

/// Observation of a landmark in an image (feature)
struct MapObservation
{
  MapObservation()
  : feat_id (UndefinedIndexT),
    frame_ptr(nullptr)
  {}
  MapObservation(IndexT f_id, Frame * frame)
  : feat_id (f_id),
    frame_ptr(frame)
  {}

  IndexT feat_id;
  Frame * frame_ptr;
};
/// Define a collection of observations of a landmark
using LandmarkObservations = Hash_Map<size_t, MapObservation>;


/// A 3D point with it's associated image observations
class MapLandmark
{
public:
  MapLandmark()
  : id_(UndefinedIndexT),
    X_(-1,-1,-1),
    normal_(0,0,0),
    ref_frame_(nullptr),
    feat_best_desc_(nullptr),
    feat_mean_scale_(0),
    active_(false),
    last_obs_step_(0),
    last_local_map_frame_id_(UndefinedIndexT)
   {}

  IndexT id_;
  Vec3 X_;
  Vec3 normal_; //mean unit vector of all viewing directions (ray that  join the point with optical center of keyframes that observe it)
  Frame * ref_frame_;

  LandmarkObservations obs_;  // map of keyframe_ids and map observation objects
  void * feat_best_desc_; // Pointer to the best descriptor of the point (most average of them all?)
  float feat_mean_scale_;  // Mean of the scales of the features associated
  size_t n_all_obs_ = 0; // Number of all observations (both frames and keyframes)

  size_t last_obs_step_;
  IndexT last_local_map_frame_id_;  // Id of frame for which the point was last added to local map
  Vec3 last_normal_;

  bool active_; // True if the point is in global map

  //Local Map  data

  // 1 - initialization point; 2- motion model/reference kf; 3-map tracking point; 4- new triangulated point
  size_t association_type_ = 0; // Through which tye of association the point was added


  const bool & isActive() const
  {
    return active_;
  }

  void setActive()
  {
    active_ = true;
  }

  const size_t getLastObsStep() const
  {
    return last_obs_step_;
  }
  void setObsStep(const size_t & step)
  {
    last_obs_step_ = step;
  }

  void updateNormal();

  const Vec3 & getWorldPosition() const
  {
    if (ref_frame_ == nullptr)
    {
      return X_;
    }
    else
    {
      // TODO: how is with relative points
      return X_;
    }
  }
  const Vec3 & getNormal() const
  {
    return normal_;
  }

  const Vec3 & getLastNormal() const
  {
    return last_normal_;
  }

  const float & getFeatureMeanScale() const
  {
    return feat_mean_scale_;
  }


  void addObservation(Frame * frame, const IndexT & feat_id);

  // Method which defines what is enough quality for the points to be added to the system
  bool isValidByConnectivityDegree(const size_t & min_degree_landmark) const;

  bool hasFrameObservation(const IndexT & frame_id);


  void setNumberOfObservations(size_t n_obs)
  {
    n_all_obs_ = n_obs;
  }
  void increaseNumberOfObservations()
  {
    n_all_obs_++;
  }
  void decreaseNumberOfObservations()
  {
    n_all_obs_ > 0 ? n_all_obs_-- : 0;
  }

};

}
}
