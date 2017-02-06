
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once


#include "openMVG/cameras/cameras.hpp"
#include "openMVG/geometry/Similarity3.hpp"
#include "openMVG/types.hpp"

using namespace openMVG;

namespace openMVG  {
namespace VSSLAM  {

class Frame;

/// Define a collection of Pose (indexed by View::id_pose)
using MapKeyFrames = Hash_Map<size_t, std::shared_ptr<Frame> >;

/// Define a collection of IntrinsicParameter (indexed by id_intrinsic)
using Intrinsics = Hash_Map<IndexT, cameras::IntrinsicBase *>;

enum MAP_POINT_TYPE
{
  EUCLIDEAN = 0,  // Absolute world reference frame
  INV_DEPTH = 1,  // Relative reference camera reference frame
  INV_DIST = 2  // Relative reference camera reference frame
};
enum MAP_CAMERA_TYPE
{
  ABSOLUTE = 0, // Absolute world reference frame
  RELATIVE = 1  // Relative reference camera reference frame
};

struct MapObservation
{
  MapObservation()
  : feat_id (std::numeric_limits<size_t>::max()),
    frame_ptr(nullptr),
    pt_ptr(nullptr),
    desc_raw_ptr(nullptr)
  {}
  MapObservation(size_t f_id, Frame * frame, Vec2 * x_ptr = nullptr, void * d_ptr = nullptr)
  : feat_id (f_id),
    frame_ptr(frame),
    pt_ptr(x_ptr),
    desc_raw_ptr(d_ptr)
  {}

  size_t feat_id;
  Frame * frame_ptr;
  Vec2 * pt_ptr;
  void * desc_raw_ptr;
};

using LandmarkObservations = Hash_Map<size_t, MapObservation>;

/// A 3D point with it's associated image observations
struct MapLandmark
{
  MapLandmark()
  : id_(std::numeric_limits<size_t>::max()),
    X_(-1,-1,-1),
    normal_(0,0,0),
    ref_frame_(nullptr),
    bestDesc_(nullptr),
    localMapFrameId_(std::numeric_limits<size_t>::max())
   {}

  size_t id_;
  Vec3 X_;
  Vec3 normal_; //mean unit vector of all viewing directions (ray that  join the point with optical center of keyframes that observe it)
  Frame * ref_frame_;

  LandmarkObservations obs_;  // map of keyframe_ids and map observation objects
  void * bestDesc_; // Pointer to the best descriptor of the point (most average of them all?)

  //Local Map  data
  size_t localMapFrameId_;  // Id of frame for which the point was last added to local map
};

using MapLandmarks = Hash_Map<size_t,MapLandmark>;

struct VSSLAM_Data
{
  // Keyframes (indexed by frame_id)
  MapKeyFrames keyframes;

  // Deque of all landmarks in the map
  MapLandmarks structure;

  // Camera intrinsics  (indexed by camera_id)
  Intrinsics cam_intrinsics;

  // Point and Camera settings
  MAP_POINT_TYPE mapPointType = MAP_POINT_TYPE::EUCLIDEAN;
  MAP_CAMERA_TYPE mapCameraType = MAP_CAMERA_TYPE::ABSOLUTE;

  size_t next_free_landmark_id = 0;

  size_t getNextFreeLandmarkId()
  {
    return next_free_landmark_id++;
  }
};



}
}
