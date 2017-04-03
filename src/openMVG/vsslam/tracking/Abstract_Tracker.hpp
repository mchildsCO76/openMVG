// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>

#include <openMVG/types.hpp>
#include <openMVG/vsslam/Frame.hpp>
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
  size_t n_keep_last_frames = 5;
  std::deque<std::shared_ptr<Frame> > mPrevFrames;  // List of

  std::shared_ptr<Frame> mLastRefFrame;
  std::shared_ptr<Frame> mPrevPrevFrame;
  std::shared_ptr<Frame> mPrevFrame;
  std::shared_ptr<Frame> mCurrentFrame;


  std::vector<std::vector<Vec2> > display_pt2d_A;
  std::vector<std::vector<Vec2> > display_pt2d_B;
  std::vector<std::vector<Vec2> > display_pt2d_C;
  std::vector<std::vector<float> > display_size_A;



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

  bool isMapInitialized()
  {
    return cartographer_->isMapInitialized();
  }

  void addPrevFrame(std::shared_ptr<Frame> & frame)
  {
    mPrevFrames.emplace_back(frame->share_ptr());
    while (mPrevFrames.size() > n_keep_last_frames)
    {
      mPrevFrames.pop_front();
    }
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
