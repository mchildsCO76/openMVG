
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once


#include "openMVG/features/features.hpp"
#include <openMVG/numeric/numeric.h>

#include <software/VSSLAM/slam/Frame.hpp>
#include <software/VSSLAM/slam/Abstract_Tracker.hpp>
#include <software/VSSLAM/slam/Abstract_FeatureExtractor.hpp>

#include <deque>
#include <set>

using namespace openMVG;
using namespace openMVG::cameras;

namespace openMVG  {
namespace VSSLAM  {

/// Store an image observation of a landmark
struct Measurement
{
  Measurement
  (
    const uint32_t & frameId,
    const Vec2f & p
  ): frameId_(frameId), pos_(p)
  { }
  Measurement( const Measurement & src ) = default ;

  uint32_t frameId_;
  Vec2f pos_;
};

/// A 3D point with it's associated image observations
struct Landmark
{
  Landmark():pt_(-1,-1,-1) {}

  Vec3 pt_;
  std::deque<Measurement> obs_;
};


/// Monocular test interface
struct SLAM_Monocular
{
  // Current frame
  std::shared_ptr<Frame> current_frame;

  // Camera
  IntrinsicBase * cam_intrinsic_;

  // Tracking
  Abstract_Tracker * tracker_;
  //

  SLAM_Monocular
  (
    Abstract_Tracker * tracker,
    IntrinsicBase * mono_cam_intrinsic
  )
  : tracker_(tracker), cam_intrinsic_(mono_cam_intrinsic)
  {
    if (tracker_)
    {
      tracker_->cam_intrinsic_ = mono_cam_intrinsic;
    }
  }

  // -------------------
  // --- System Initialization
  // -------------------
  bool isReady()
  {
    // if either of tracker or intrinsic are not initialized is not ready
    if (!tracker_)
    {
      std::cerr << "ERROR: MonoSLAM: Tracker not initialized!" << std::endl;
      return false;
    }
    if (!cam_intrinsic_)
    {
      std::cerr << "ERROR: MonoSLAM: Camera intrinsics not initialized!" << std::endl;
      return false;
    }
    return true;
  }



  bool nextFrame
  (
    const image::Image<unsigned char> & ima,
    const size_t frameId
  )
  {
    std::cout<<"Frame "<<frameId<<"\n";
    // Create Frame
    current_frame = std::make_shared<Frame>(frameId);

    double startTime = omp_get_wtime();
    // Track frame
    tracker_->track(ima,current_frame);

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Track time:"<<secsElapsed<<"\n";

    tracker_->printTrackingStatus();

  }




};

}
}
