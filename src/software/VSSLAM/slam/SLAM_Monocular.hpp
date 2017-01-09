
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

  // Tracking
  Abstract_Tracker * tracker_;
  //


  SLAM_Monocular
  (
    Abstract_Tracker * tracker,
    const uint32_t maxTrackedFeatures = 1500
    // Add an abstract camera model here
  )
  : tracker_(tracker)
  {
    tracker_->max_tracked_points = maxTrackedFeatures;
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

    const Abstract_Tracker::TRACKING_STATUS systemStatus = tracker_->trackingStatus;
    switch (systemStatus)
    {
    case Abstract_Tracker::TRACKING_STATUS::OK:
      std::cout<<"Tracking OK\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::LOST:
      std::cout<<"Tracking LOST\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::NOT_INIT:
      std::cout<<"Tracking NOT_INIT\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::INIT:
      std::cout<<"Tracking INIT\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::IDLE:
      std::cout<<"Tracking IDLE\n";
      break;
    }

    /*
    // Check if system is initialized
    if (systemStatus == Abstract_Tracker::TRACKING_STATUS::INIT)
    {
      MonocularSystemInitialization(ima);
    }*/
  }

/*
  // Try to initialize system
  bool MonocularSystemInitialization
  (
    const image::Image<unsigned char> & ima
  )
  {
    std::cout<<"System not initialized yet!\n";

    // Check if we have enough features in the frame
    if (tracker_->trackingStatus != Abstract_Tracker::TRACKING_STATUS::FAILED)
    {
      // Check if we have a reference init frame
      if (!init_ref_frame)
      {
        // Set current frame as the new reference frame for initialization
        std::cout<<"Set: INIT REFERENCE\n";
        init_ref_frame.reset(current_frame.get());
      }
      else if (tracker_->trackingStatus == Abstract_Tracker::TRACKING_STATUS::OK)
      {
        // Tracking from previous frame was successful
        std::cout<<"Try to initialize\n";

      }
    }
  }
*/


};

}
}
