// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>

#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/types.hpp>
#include <openMVG/vsslam/mapping/Cartographer.hpp>
#include <openMVG/vsslam/tracking/MotionModel.hpp>

using namespace openMVG;
using namespace openMVG::cameras;

namespace openMVG  {
namespace VSSLAM  {

class Abstract_Tracker
{
public:
  enum TRACKING_STATUS
  {
    NOT_INIT = 0,
    INIT = 1,
    OK = 2,
    LOST = 3,
    IDLE = 4
  };
protected:
  // Map
  Cartographer * cartographer_;

  // Tracking stats
  TRACKING_STATUS trackingStatus = TRACKING_STATUS::NOT_INIT;

  // Motion model
  MotionModel motionModel;

  // Initialization
  std::shared_ptr<Frame> init_ref_frame;


public:

  // ---------------
  // Parameters
  // ---------------

  // Frames
  std::shared_ptr<Frame> mLastRefFrame;
  std::shared_ptr<Frame> mPrevFrame;
  std::shared_ptr<Frame> mCurrentFrame;



  // ---------------
  // Methods
  // ---------------

  virtual ~Abstract_Tracker(){};

  // Try to track the pt_to_track features
  virtual bool track
  (
    const image::Image<unsigned char> & ima,
    std::shared_ptr<Frame> current_frame
  ) = 0;

  void setCartographer(Cartographer * cart)
  {
    cartographer_ = cart;
  }

  TRACKING_STATUS getTrackingStatus()
  {
    return trackingStatus;
  }

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
