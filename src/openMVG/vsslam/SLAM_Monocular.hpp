
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once


#include <openMVG/types.hpp>
#include "openMVG/features/features.hpp"
#include <openMVG/numeric/numeric.h>
#include <openMVG/vsslam/Camera.hpp>
#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/tracking/Abstract_Tracker.hpp>
#include <openMVG/vsslam/tracking/Abstract_FeatureExtractor.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include <openMVG/vsslam/mapping/Cartographer.hpp>
#include <deque>
#include <set>

#include <iostream>
#include <memory>

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
  Hash_Map<size_t, std::shared_ptr<Camera> > cameras;
  //IntrinsicBase * cam_intrinsic_;

  // Tracking
  Abstract_Tracker * tracker_;

  // Map
  std::shared_ptr<Cartographer> cartographer_;

  SLAM_Monocular
  (
    Abstract_Tracker * tracker
  )
  : tracker_(tracker)
  {
    cartographer_ = std::make_shared<Cartographer>();
    if (tracker_)
    {
      tracker_->cartographer_ = cartographer_.get();
    }
  }

  void setFeatureExtractor(Abstract_FeatureExtractor * f_extractor)
  {
    cartographer_->feature_extractor_ = f_extractor;
  }

  int createCamera(const CameraParams & cam_params)
  {
    std::shared_ptr<Camera> cam = std::make_shared<Camera>();
    cam->bCalibrated = cam_params.bCalibrated;
    switch(cam_params.camera_model)
    {
      case PINHOLE_CAMERA:
        cam->cam_intrinsic_ = std::make_shared<Pinhole_Intrinsic>
          (cam_params.img_width, cam_params.img_height, cam_params.focal, cam_params.ppx, cam_params.ppy);
        cam->cam_intrinsic_ptr = cam->cam_intrinsic_.get();
      break;
      case PINHOLE_CAMERA_RADIAL1:
        cam->cam_intrinsic_ = std::make_shared<Pinhole_Intrinsic_Radial_K1>
          (cam_params.img_width, cam_params.img_height, cam_params.focal, cam_params.ppx, cam_params.ppy, 0.0); // setup no distortion as initial guess
      break;
      case PINHOLE_CAMERA_RADIAL3:
        cam->cam_intrinsic_ = std::make_shared<Pinhole_Intrinsic_Radial_K3>
          (cam_params.img_width, cam_params.img_height, cam_params.focal, cam_params.ppx, cam_params.ppy, 0.0, 0.0, 0.0);  // setup no distortion as initial guess
      break;
      case PINHOLE_CAMERA_BROWN:
        cam->cam_intrinsic_ =std::make_shared<Pinhole_Intrinsic_Brown_T2>
          (cam_params.img_width, cam_params.img_height, cam_params.focal, cam_params.ppx, cam_params.ppy, 0.0, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
      break;
      case PINHOLE_CAMERA_FISHEYE:
        cam->cam_intrinsic_ =std::make_shared<Pinhole_Intrinsic_Fisheye>
          (cam_params.img_width, cam_params.img_height, cam_params.focal, cam_params.ppx, cam_params.ppy, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
      break;
      default:
        return -1;
    }
    if (cam_params.camera_model != PINHOLE_CAMERA && cam_params.bCalibrated)
    {
      cam->cam_intrinsic_undist = std::make_shared<Pinhole_Intrinsic>
      (cam_params.img_width, cam_params.img_height, cam_params.focal, cam_params.ppx, cam_params.ppy);
      cam->cam_intrinsic_ptr = cam->cam_intrinsic_undist.get();
    }
    else
    {
      cam->cam_intrinsic_ptr = cam->cam_intrinsic_.get();
    }
    // Insert camera into database
    cameras[cameras.size()]=cam;
    return cameras.size();
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
    if (cameras.size()<1)
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
    const size_t frameId,
    const size_t camId
  )
  {
    std::cout<<"Frame "<<frameId<<"\n";
    // Create Frame
    current_frame = std::make_shared<Frame>(frameId, camId, cameras[camId].get());

    double startTime = omp_get_wtime();
    // Track frame
    bool bTrack = tracker_->track(ima,current_frame);

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Track time:"<<secsElapsed<<"\n";

    if (!bTrack)
    {
      std::cout<<"TRY RELOCALIZATION!\n";
    }




    tracker_->printTrackingStatus();

  }




};

}
}
