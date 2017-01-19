
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once


#include "openMVG/features/features.hpp"
#include <openMVG/numeric/numeric.h>

#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/tracking/Abstract_Tracker.hpp>
#include <openMVG/vsslam/tracking/Abstract_FeatureExtractor.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include <openMVG/vsslam/mapping/Cartographer.hpp>
#include <deque>
#include <set>

using namespace openMVG;
using namespace openMVG::cameras;

namespace openMVG  {
namespace VSSLAM  {



/// Monocular test interface
struct SLAM_Monocular
{
  // Current frame
  std::shared_ptr<Frame> current_frame;

  // Camera
  IntrinsicBase * cam_intrinsic_;

  // Tracking
  Abstract_Tracker * tracker_;

  // Map
  std::shared_ptr<Cartographer> cartographer_;

  SLAM_Monocular
  (
    Abstract_Tracker * tracker,
    IntrinsicBase * mono_cam_intrinsic
  )
  : tracker_(tracker), cam_intrinsic_(mono_cam_intrinsic)
  {
    cartographer_ = std::make_shared<Cartographer>();
    if (tracker_)
    {
      tracker_->cam_intrinsic_ = mono_cam_intrinsic;
      tracker_->cartographer_ = cartographer_.get();
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

  void addCameraIntrinsics(const size_t & cam_id, IntrinsicBase * cam_intrinsic_)
  {
    cartographer_->addCameraIntrinsicData(cam_id,cam_intrinsic_);
  }


  bool nextFrame
  (
    const image::Image<unsigned char> & ima,
    const size_t frameId
  )
  {
    std::cout<<"Frame "<<frameId<<"\n";
    // Create Frame
    current_frame = std::make_shared<Frame>(frameId, 0, cam_intrinsic_);

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
