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


    if (time_data.b_export_stats_file)
    {
      time_data.stats_file.open(time_data.s_stats_file_path.c_str(),std::ios::out );
      time_data.stats_file << "frame_id;"
          <<" keyframe;"
          <<" keyframe_reason;"
          <<" track_total;"
          <<" track_detect;"
          <<" track_track;"

          <<" track_mm;"
          <<" track_mm_poseBA;"
          <<" track_rf;"
          <<" track_rf_poseBA;"
          <<" track_rm;"
          <<" track_rm_poseBA;"
          <<" track_lm;"
          <<" track_lm_poseBA;"

          <<" track_new_pts;"
          <<" track_new_pts_localBA;"
          <<" track_global_BA;"

          <<" matches_mm;"
          <<" matches_mm_outliers;"
          <<" matches_rf;"
          <<" matches_rf_outliers;"
          <<" matches_rm;"
          <<" matches_lm;"

          <<" map_global_frames;"
          <<" map_global_landmarks;"
          <<" map_local_landmarks;"

          <<" map_added_global_landmarks;"
          <<" map_added_local_landmarks;"
          <<" map_added_local_to_global_landmarks;"
          <<" map_removed_local_landmarks_outliers;"
          <<" map_removed_local_landmarks_inactive;"
          <<"\n";
      time_data.stats_file.flush();

    }

  }
  SLAM_System::~SLAM_System()
  {

    if (time_data.b_export_stats_file)
    {
      time_data.stats_file.close();
    }
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

    if (time_data.b_enable_time_stats || time_data.b_enable_features_stats)
    {
      time_data.restartData();
      time_data.frame_id = id_frame;
    }

    // Create frame
    Camera * ptr_cam = map_cameras_[id_cam].get();
    frame_current_ = std::make_shared<Frame>(id_frame, time_frame, ptr_cam);


    // Time statistics
    time_data.startTimer(time_data.d_track_track);

    // Track frame
    tracker_->track(ima,frame_current_,ptr_cam->getMaskImagePtr());

    cartographer_->setMapStats(time_data);

    // Time statistics
    time_data.stopTimer(time_data.d_track_track);

    // Show tracking status
    tracker_->printTrackingStatus();

    if (time_data.b_enable_time_stats || time_data.b_enable_features_stats)
    {
      exportStatistics(time_data);
      printStatistics(time_data);
    }
  }

  void SLAM_System::exportStatistics(VSSLAM_Time_Stats & stats)
  {
    if (time_data.b_export_stats_file)
    {
      time_data.stats_file
         << stats.frame_id<<"; "
         << stats.b_keyframe << "; "
         << stats.keyframe_reason << "; "

         << stats.d_track_track << "; "
         << stats.d_feat_detection << "; "
         << stats.d_feat_tracking << "; "

         << stats.d_feat_tracking_mm << "; "
         << stats.d_feat_pose_opt_mm << "; "
         << stats.d_feat_tracking_rf << "; "
         << stats.d_feat_pose_opt_rf << "; "
         << stats.d_feat_tracking_rm << "; "
         << stats.d_feat_pose_opt_rm << "; "
         << stats.d_feat_tracking_lm << "; "
         << stats.d_feat_pose_opt_lm << "; "

         << stats.d_feat_new_pts << "; "
         << stats.d_feat_pose_local << "; "
         << stats.d_feat_add_to_global << "; "

         << stats.i_matches_mm << "; "
         << stats.i_matches_mm_outliers << "; "
         << stats.i_matches_rf << "; "
         << stats.i_matches_rf_outliers << "; "
         << stats.i_matches_rm << "; "
         << stats.i_matches_lm << "; "

         << stats.global_frames << "; "
         << stats.global_landmarks << "; "
         << stats.local_landmarks << "; "

         << stats.added_global_landmarks << "; "
         << stats.added_local_landmarks << "; "
         << stats.added_local_to_global_landmarks << "; "
         << stats.removed_local_landmarks_outliers << "; "
         << stats.removed_local_landmarks_inactive << "; "

         <<"\n";
     time_data.stats_file.flush();

    }
  }


  void SLAM_System::printStatistics(VSSLAM_Time_Stats & stats)
  {
    std::cout<<"---------------------------\n"
        <<"----------STATS-------------\n"
        <<"---------------------------\n";
    std::cout<<"Frame ID: "<<stats.frame_id<<"\n"
        <<"Keyframe: "<<stats.b_keyframe<<"\n"
        <<"Keyframe reason: "<<stats.keyframe_reason<<"\n"


        <<"--------Times-----------\n"
        <<"Tracking Total: "<<stats.d_track_track<<"\n"
        <<"Detection: "<<stats.d_feat_detection<<"\n"
        <<"Tracking: "<<stats.d_feat_tracking<<"\n"
        <<"Tracking MM: "<<stats.d_feat_tracking_mm<<"\n"
        <<"Pose MM: "<<stats.d_feat_pose_opt_mm<<"\n"

        <<"Tracking RF: "<<stats.d_feat_tracking_rf<<"\n"
        <<"Pose RF: "<<stats.d_feat_pose_opt_rf<<"\n"
        <<"Tracking RM: "<<stats.d_feat_tracking_rm<<"\n"
        <<"Pose RM: "<<stats.d_feat_pose_opt_rm<<"\n"
        <<"Tracking LM: "<<stats.d_feat_tracking_lm<<"\n"
        <<"Pose LM: "<<stats.d_feat_pose_opt_lm<<"\n"
        <<"Tracking N_PTS: "<<stats.d_feat_new_pts<<"\n"
        <<"Pose Local: "<<stats.d_feat_pose_local<<"\n"
        <<"Global BA: "<<stats.d_feat_add_to_global<<"\n"

        <<"--------Features-----------\n"
        <<"Matches MM / Outliers: "<<stats.i_matches_mm<<"/"<<stats.i_matches_mm_outliers<<"\n"
        <<"Matches RF / Outliers: "<<stats.i_matches_rf<<"/"<<stats.i_matches_rf_outliers<<"\n"
        <<"Matches RM: "<<stats.i_matches_rm<<"\n"
        <<"Matches LM: "<<stats.i_matches_lm<<"\n"

        <<"--------Map-----------\n"
        <<"Global Frames: "<<stats.global_frames<<"\n"
        <<"Global Landmarks: "<<stats.global_landmarks<<"\n"
        <<"Local Landmarks: "<<stats.local_landmarks<<"\n"
        <<"-------------------\n"
        <<"Added Global Landmarks: "<<stats.added_global_landmarks<<"\n"
        <<"Added Local Landmarks: "<<stats.added_local_landmarks<<"\n"
        <<"Added Local-Global Landmarks: "<<stats.added_local_to_global_landmarks<<"\n"
        <<"Removed Local Landmarks (outliers): "<<stats.removed_local_landmarks_outliers<<"\n"
        <<"Removed Local Landmarks (inactive): "<<stats.removed_local_landmarks_inactive<<"\n"
        <<"-------------------\n";




  }

}
}
