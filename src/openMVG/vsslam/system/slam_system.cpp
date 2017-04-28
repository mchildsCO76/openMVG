// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <memory>

#include <openMVG/vsslam/system/slam_system.hpp>
#include <openMVG/vsslam/tracking/Tracker_Features.hpp>

namespace openMVG {
namespace vsslam {
  // Construction
  SLAM_System::SLAM_System(std::shared_ptr<VSSLAM_Parameters> & params)
  {
    std::cout<<"VSSLAM: [System] Created new object\n";
    params_ = params->share_ptr();
  }
  void SLAM_System::setTracker(std::unique_ptr<Abstract_Tracker> & tracker)
  {
    tracker_ = std::move(tracker);
    if (feat_extractor_)
    {
      dynamic_cast<Tracker_Features *>(tracker_.get())->setFeatureExtractor(feat_extractor_.get());
      std::cout<<"VSSLAM: [System] Tracker assigned - feat extractor\n";
    }
    if (feat_matcher_)
    {
      dynamic_cast<Tracker_Features *>(tracker_.get())->setFeatureMatcher(feat_matcher_.get());
      std::cout<<"VSSLAM: [System] Tracker assigned - feat matcher\n";
    }
    if (cartographer_)
    {
      tracker_->setCartographer(cartographer_.get());
      std::cout<<"VSSLAM: [System] Tracker assigned - cartographer\n";
    }
    std::cout<<"VSSLAM: [System] Tracker assigned\n";
  }

  void SLAM_System::setFeatureExtractor(std::unique_ptr<Abstract_Feature_Extractor> & extractor)
  {
    feat_extractor_ = std::move(extractor);
    if (tracker_)
      dynamic_cast<Tracker_Features *>(tracker_.get())->setFeatureExtractor(feat_extractor_.get());
    if (cartographer_)
      cartographer_->setFeatureExtractor(feat_extractor_.get());
    std::cout<<"VSSLAM: [System] Feature extractor assigned\n";
  }
  void SLAM_System::setFeatureMatcher(std::unique_ptr<Abstract_Feature_Matcher> & matcher)
  {
    feat_matcher_ = std::move(matcher);
    if (tracker_)
      dynamic_cast<Tracker_Features *>(tracker_.get())->setFeatureMatcher(feat_matcher_.get());
    std::cout<<"VSSLAM: [System] Feature matcher assigned\n";
  }

  bool SLAM_System::createCartographer
  (
    MAP_FRAME_TYPE map_frame_type,
    MAP_LANDMARK_TYPE map_landmark_type,
    MAP_OPTIMIZATION_TYPE global_BA_type,
    MAP_OPTIMIZATION_TYPE local_BA_type
  )
  {
    cartographer_.reset(new Cartographer(params_,map_frame_type,map_landmark_type,global_BA_type,local_BA_type));

    if (feat_extractor_)
      cartographer_->setFeatureExtractor(feat_extractor_.get());

    if (tracker_)
    {
      tracker_->setCartographer(cartographer_.get());
    }
    std::cout<<"VSSLAM: [System] Cartographer initialization OK\n";
    return true;
  }

  bool SLAM_System::isReady()
  {
    // if either of tracker or intrinsic are not initialized is not ready
    if (!tracker_ || !tracker_->isReady())
    {
      std::cerr << "VSSLAM: [System] Tracker is not assigned\n";
      return false;
    }
    if (map_cameras_.empty())
    {
      std::cerr << "VSSLAM: [System] No cameras are initialized\n";
      return false;
    }
    if (!cartographer_ || !cartographer_->isReady())
    {
      std::cerr << "VSSLAM: [System] Cartographer is not assigned\n";
      return false;
    }
    return true;
  }

  // Insert a camera together with a mask image for processing
  IndexT SLAM_System::createCamera(const CameraParameters & param_cam)
  {
    // Define new id for camera
    IndexT id_new_cam = map_cameras_.size();
    // Create new camera
    std::shared_ptr<Camera> cam = std::make_shared<Camera>(id_new_cam, param_cam.b_calibrated);
    // Create intrinsic parameters
    switch ( param_cam.camera_model)
    {
      case PINHOLE_CAMERA:
        cam->ptr_intrinsic_ = std::make_shared<Pinhole_Intrinsic>
        (param_cam.img_width, param_cam.img_height, param_cam.focal, param_cam.ppx, param_cam.ppy);
      break;
      case PINHOLE_CAMERA_RADIAL1:
        cam->ptr_intrinsic_ = std::make_shared<Pinhole_Intrinsic_Radial_K1>
        (param_cam.img_width, param_cam.img_height, param_cam.focal, param_cam.ppx, param_cam.ppy, 0.0); // setup no distortion as initial guess
      break;
      case PINHOLE_CAMERA_RADIAL3:
        cam->ptr_intrinsic_ = std::make_shared<Pinhole_Intrinsic_Radial_K3>
              (param_cam.img_width, param_cam.img_height, param_cam.focal, param_cam.ppx, param_cam.ppy, 0.0, 0.0, 0.0); // setup no distortion as initial guess
      break;
      case PINHOLE_CAMERA_BROWN:
        cam->ptr_intrinsic_ = std::make_shared<Pinhole_Intrinsic_Brown_T2>
              (param_cam.img_width, param_cam.img_height, param_cam.focal, param_cam.ppx, param_cam.ppy, 0.0, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
      break;
      case PINHOLE_CAMERA_FISHEYE:
        cam->ptr_intrinsic_ = std::make_shared<Pinhole_Intrinsic_Fisheye>
              (param_cam.img_width, param_cam.img_height, param_cam.focal, param_cam.ppx, param_cam.ppy, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
      break;
      default:
        return std::numeric_limits<IndexT>::max();
    }

    // If camera has distortion is calibrated create an undistorted intriniscs
    if (param_cam.camera_model != PINHOLE_CAMERA && param_cam.b_calibrated)
    {
      cam->ptr_intrinsic_undist_ = std::make_shared<Pinhole_Intrinsic>
      (param_cam.img_width, param_cam.img_height, param_cam.focal, param_cam.ppx, param_cam.ppy);
      cam->ptr_intrinsic_valid_ = cam->ptr_intrinsic_undist_.get();
    }
    else
    {
      cam->ptr_intrinsic_valid_ = cam->ptr_intrinsic_.get();
    }

    // Compute camera borders
    cam->computeImageBorders();

    // Insert camera into database
    map_cameras_[id_new_cam] = cam;

    std::cout<<"VSSLAM: [System] Camera "<<id_new_cam<<" created\n";
    return id_new_cam;
  }

  // Insert a camera together with a mask image for processing
  IndexT SLAM_System::createCamera(const CameraParameters & param_cam, image::Image<unsigned char> & mask /*= nullptr*/)
  {
    // Create new camera
    IndexT id_cam = createCamera(param_cam);

    // Check if camera is valid
    if (id_cam == std::numeric_limits<IndexT>::max())
      return id_cam;

    // Save mask image
    map_cameras_[id_cam]->setMaskImage(mask);

    std::cout<<"VSSLAM [System] Camera "<<id_cam<<" added mask image\n";

    return id_cam;
  }

  bool SLAM_System::addMaskImageToCamera(const IndexT & id_cam, image::Image<unsigned char> & mask)
  {
    // Wrong id of camera
    if (map_cameras_.find(id_cam) == map_cameras_.end())
      return false;
    else
    {
      map_cameras_[id_cam]->setMaskImage(mask);
      std::cout<<"VSSLAM [System] Camera "<<id_cam<<" added mask image\n";
    }
    return true;
  }


  // Process images
  bool SLAM_System::nextFrame
  (
    const image::Image<unsigned char> & ima,
    const IndexT & id_frame,
    const IndexT & id_cam,
    const double & time_frame
  )
  {
    std::cout<<"VSSLAM [System] Frame "<<id_frame<<" of camera: "<< id_cam<<" at time: "<<time_frame<<"\n";

    // Create frame
    Camera * ptr_cam = map_cameras_[id_cam].get();
    frame_current_ = std::make_shared<Frame>(id_frame, time_frame, ptr_cam);

    // Track frame
    tracker_->track(ima,frame_current_,ptr_cam->getMaskImagePtr());

    // Show tracking status
    tracker_->printTrackingStatus();

  }


}
}
