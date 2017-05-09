// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include <openMVG/types.hpp>
#include "openMVG/stl/split.hpp"
//#include "openMVG/image/image.hpp" // swine
#include <openMVG/cameras/cameras.hpp>


#include <openMVG/vsslam/vsslam_parameters.hpp>
#include <openMVG/vsslam/system/Frame.hpp>
#include <openMVG/vsslam/vsslam_data.hpp>
#include <openMVG/vsslam/mapping/MapLandmark.hpp>

namespace openMVG {
namespace vsslam {

struct VSSLAM_Map
{
  MAP_LANDMARK_TYPE map_landmark_type_ = MAP_LANDMARK_TYPE::GLOBAL_EUCLIDEAN;
  MAP_FRAME_TYPE map_frame_type_ = MAP_FRAME_TYPE::GLOBAL;

  // Keyframes (indexed by frame_id)
  MapFrames map_frame_;
  // Landmarks (indexed by landmark_id)
  MapLandmarks map_landmarks_;
  // Camera intrinsics  (indexed by camera_id)
  Intrinsics map_intrinsics;

  // Keep track of used landmarkIds
  IndexT next_free_landmark_id = 0;

  void clear()
  {
    map_frame_.clear();
    map_landmarks_.clear();
    map_intrinsics.clear();
    next_free_landmark_id = 0;
  }

  size_t getNumberOfFrames()
  {
    return map_frame_.size();
  }

  size_t getNumberOfLandmarks()
  {
    return map_landmarks_.size();
  }

  IndexT getNextFreeLandmarkId()
  {
    return next_free_landmark_id++;
  }
  bool addFrame(const std::shared_ptr<Frame> & frame)
  {
    // Return false if frame already exists in the map
    if (map_frame_.find(frame->getFrameId()) != map_frame_.end())
      return false;

    map_frame_[frame->getFrameId()] = frame->share_ptr();

    if (map_intrinsics.find(frame->getCamId()) != map_intrinsics.end())
    {
      map_intrinsics[frame->getCamId()] = frame->getCameraIntrinsics();
    }

    return true;
  }

  MapLandmark * addLandmark(std::unique_ptr<MapLandmark> & landmark_new)
  {
    // Get new ID for landmark
    const IndexT id_landmark = getNextFreeLandmarkId();
    landmark_new->id_ = id_landmark;

    // Add landmark to structure
    MapLandmark * ptr_map_landmark = landmark_new.get();
    map_landmarks_[id_landmark] = std::move(landmark_new);

    return ptr_map_landmark;
  }

  std::unique_ptr<MapLandmark> & getLandmark(const IndexT & id_landmark)
  {
    return map_landmarks_[id_landmark];
  }

  void removeLandmark(IndexT & id_landmark)
  {
    map_landmarks_.erase(id_landmark);
  }

  void removeFrame(const IndexT & id_frame)
  {
    map_frame_.erase(id_frame);
  }

};

}
}
