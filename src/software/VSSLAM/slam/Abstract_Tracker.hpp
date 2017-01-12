// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>

#include <software/VSSLAM/slam/Frame.hpp>
#include <openMVG/types.hpp>

using namespace openMVG;
using namespace openMVG::cameras;

namespace openMVG  {
namespace VSSLAM  {


struct Abstract_Tracker
{
  enum TRACKING_STATUS
  {
    NOT_INIT = 0,
    INIT = 1,
    OK = 2,
    LOST = 3,
    IDLE = 4
  };

  // Camera intrinsics
  IntrinsicBase * cam_intrinsic_;

  // Frames
  std::shared_ptr<Frame> mLastRefFrame;
  std::shared_ptr<Frame> mPrevFrame;
  std::shared_ptr<Frame> mCurrentFrame;

  // Initialization
  std::shared_ptr<Frame> init_ref_frame;

  // Tracking stats
  TRACKING_STATUS trackingStatus = TRACKING_STATUS::NOT_INIT;

  Abstract_Tracker() = default;

  // Try to track the pt_to_track features
  virtual bool track
  (
    const image::Image<unsigned char> & ima,
    std::shared_ptr<Frame> current_frame
  ) = 0;

  // suggest new feature point for tracking (count point are kept)
  /*virtual bool detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count,
    const size_t min_count = 0
  ) const = 0;*/

  void printTrackingStatus()
  {
    std::cout<<"Tracking STATUS: ";
    switch (trackingStatus)
    {
    case Abstract_Tracker::TRACKING_STATUS::OK:
      std::cout<<"OK\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::LOST:
      std::cout<<"LOST\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::INIT:
      std::cout<<"INIT\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::NOT_INIT:
      std::cout<<"NOT_INIT\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::IDLE:
      std::cout<<"IDLE\n";
      break;
    }
  }
};

} // namespace VO
} // namespace openMVG
