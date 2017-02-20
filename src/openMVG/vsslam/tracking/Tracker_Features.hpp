// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>
#include <tuple>

#include <openMVG/types.hpp>
#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/tracking/Abstract_Tracker.hpp>
#include <openMVG/vsslam/tracking/Abstract_FeatureExtractor.hpp>

#include <openMVG/vsslam/tracking/PoseEstimation.hpp>
#include <openMVG/vsslam/tracking/Matching.hpp>
#include "openMVG/multiview/triangulation.hpp"

namespace openMVG  {
namespace VSSLAM  {

class Tracker_Features : public Abstract_Tracker
{
private:
  /// Feature extractor
  Abstract_FeatureExtractor * featureExtractor_;

  // ---------------
  // Parameters
  // ---------------
  size_t max_tracked_points = 1500;


  // ---------------
  // Tracking init settings
  // ---------------
  size_t init_track_min_matches = 30; // Min number of points needed for initialization of tracking
  float init_min_cos_angle_pt = 0.99995; // Min angle between rays for the point to be triangulated (0.99995 ~ 0.5deg; 0.99998 ~ 0.36deg; 0.9998 ~ 1.15deg)

  float default_reproj_thresh_2 = 16.0; // squared reprojection error used as default value for frames
  // ---------------
  // Matching settings
  // ---------------
  float init_match_desc_ratio_2 = 0.8*0.8; // Matching ration of descriptors
  float track_match_desc_ratio_2 = 0.8*0.8;

  size_t track_min_matches = 10;  // Min matches needed for the tracking to be successful
  float track_mm_win_size = 20*20;

  size_t track_local_map_size = 5;  // Number of best connected frames are used for constructing a local network

  //size_t init_min_map_points = 50; // Minimum points needed for map initialization
  //size_t init_min_pts = 50; // Min points needed for successful initialization
// Min points needed for tracking between frames

  float track_epipolar_desc_ratio = 0.6*0.6;

  float track_local_map_win_size = 20*20;

  size_t triangule_local_map_size = 5;


  // ---------------
  // Tracking settings
  // ---------------
  size_t min_init_ref_tracks = 500; // Min number of tracks reference frame for initialization has to have
  size_t min_frame_tracks = 300; // Min number of tracks frame has to have
  size_t max_frame_tracks = 0; // Max number of feats detected in a frame (0 - unlimited)



public:
  // ---------------
  // Parameters
  // ---------------
  Tracker_Features
  (
    Abstract_FeatureExtractor * featExtractor,
    const size_t max_features_tracked = 1500
  ): featureExtractor_(featExtractor), max_tracked_points(max_features_tracked)
  {}

  ~Tracker_Features() = default;

  void setFeatureExtractor(Abstract_FeatureExtractor * featExtractor)
  {
    featureExtractor_ = featExtractor;
  }

  void setMaxFeaturesTracked(size_t max_feats)
  {
    max_tracked_points = max_feats;
  }


  /// Try to track current point set in the provided image
  /// return false when tracking failed (=> to send frame to relocalization)
  bool track
  (
    const image::Image<unsigned char> & ima,
    std::shared_ptr<Frame> current_frame
  ) override
  {
    bool track_status = false;
    // Set current frame
    mCurrentFrame = current_frame->share_ptr();



    // -------------------------
    // -- Detect features
    // -------------------------

    std::cout<<"Tracker: Start feature detection\n";
    double startTime = omp_get_wtime();

    if (!detect(ima,mCurrentFrame.get(),max_tracked_points,0))
    {
      // No features detected on the frame
      std::cout<<"Tracker: Error - no features detected!\n";
      return false;
    }
    const size_t & n_feats_detected = mCurrentFrame->getNumberOfFeatures();
    std::cout<<"Tracker: Features detected: "<<n_feats_detected<<"\n";

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Tracker: End feature detection ("<<secsElapsed<<" s)\n";


    // -------------------------
    // -- Tracking
    // -------------------------

    // If tracking is not initialized
    if (trackingStatus == TRACKING_STATUS::NOT_INIT)
    {
      std::cout << "Tracker: Try system initialization! Frame id: " << mCurrentFrame->getFrameId()<<"\n";
      // Check if enough features are detected
      if (n_feats_detected > init_track_min_matches)
      {
        track_status = trySystemInitialization(ima);
      }
      else
      {
        // not enough features detected
        std::cout<<"Tracker: Insufficient number of features detected!\n";
        resetSystemInitialization();
        trackingStatus = Abstract_Tracker::TRACKING_STATUS::NOT_INIT;
        track_status = false;
      }
    }
    else if (trackingStatus == TRACKING_STATUS::INIT)
    {
      std::cout << "Tracker: Try system initialization! Frame id: " << mCurrentFrame->getFrameId()<<" with Ref frame id: "<<init_ref_frame->getFrameId()<<"\n";
      // Check if enough features are detected
      if (n_feats_detected > init_track_min_matches)
      {
        track_status = trySystemInitialization(ima);
      }
      else
      {
        // not enough features detected - we try with next frame (keeping init reference frame)
        std::cout<<"Tracker: Insufficient number of features detected!\n";
        track_status = false;
      }
    }
    else if (trackingStatus == TRACKING_STATUS::OK)
    {
      bool b_track_OK = false;

      // ----------------
      // -- Motion Model tracking
      // ----------------
      // Use previous tracking to predict where the pose will be and try matching with previously triangulated map points
      if (motionModel.isValid())
      {
        std::cout<<"Tracker: Valid Motion model!\n";

        // Match current frame with lastFrame
        b_track_OK = trackWithMotionModel();
      }
      std::cout<<"Tracker: Tracking with MM: "<<b_track_OK<<"\n";

      // If motion model tracking didnt work we try with reference frame
      // ----------------
      // -- Reference frame tracking
      // ----------------
      if (!b_track_OK)
      {
        b_track_OK = trackWithReferenceFrame();
        std::cout<<"Tracker: Tracking with RF: "<<b_track_OK<<"\n";
      }

      std::vector<std::unique_ptr<MapLandmark> > vec_new_pts_3D;

      if (b_track_OK)
      {
        // ----------------
        // -- Local map tracking
        // ----------------
        trackLocalMap();

        // ----------------
        // -- Find new points that we can triangulate
        // ----------------
        std::vector<Frame *> local_map_frames;
        triangulateNewPoints(local_map_frames, vec_new_pts_3D);

        for (auto & pt : vec_new_pts_3D)
        {
          if (pt->obs_.size()==0)
            std::cout<<"KKKKKKKKKKKKKKKKKKK\n";
        }
        // ----------------
        // -- Perform BA on all the matches existing points and new trinagulated points
        // ----------------
        VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
        options.linear_solver_type_ = ceres::DENSE_SCHUR;
        VSSLAM_Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        local_map_frames.push_back(mCurrentFrame.get());
        if (!bundle_adjustment_obj.OptimizePose(&local_map_frames, mCurrentFrame.get(), nullptr, &vec_new_pts_3D))
        {
          return false;
        }
        for (auto & pt : vec_new_pts_3D)
        {
          if (pt->obs_.size()==0)
            std::cout<<"VVVVVVVVVVVVVVV\n";
        }
        removeOutliersFromFrames(local_map_frames,vec_new_pts_3D);
        size_t n_pts_frame = mCurrentFrame->getNumberMapPoints();

        for (auto & pt : vec_new_pts_3D)
        {
          if (pt->obs_.size()==0)
            std::cout<<"LLLLLLLLLLLLLLLLLLLLL\n";
        }


        std::cout<<"Tracker: New Trianguulated points: "<<vec_new_pts_3D.size()<<"\n";

        if (!cartographer_->isMapInitialized())
        {
          if (n_pts_frame < init_track_min_matches)
          {
            b_track_OK = false;
          }
        }
        else
        {
          if (n_pts_frame < track_min_matches)
          {
            b_track_OK = false;
          }
        }
      }

      if (b_track_OK)
      {
        // ----------------
        // -- Decide if current frame is keyframe and if should be added to the system
        // ----------------
        bool bKeyframe = true;

        if (bKeyframe)
        {
          if (!cartographer_->isMapInitialized())
          {
            // Enough points to add to initialization map
            if(!cartographer_->initializationAddStep(mCurrentFrame, &vec_new_pts_3D))
            {
              // What happens if initialization of map fails
              // reset tracker
            }
          }
          else
          {
            cartographer_->addStep(mCurrentFrame, &vec_new_pts_3D);
          }
          mCurrentFrame->computeFrameVisibilityConnections();
          mLastRefFrame = mCurrentFrame->share_ptr();
        }

        // Set motion model
        motionModel.updateMotionModel(mPrevFrame.get(),mCurrentFrame.get());

        track_status = true;
      }
      else
      {
        /*
        // We couldnt match with previous frame
        if (cartographer_->getNumberOfKeyframes() < 5)
        {
          std::cout<<"Tracker: System lost in intial 5 frames! New init required\n";
          // Set tracking as not initialized
          resetSystemInitialization();
          trackingStatus = Abstract_Tracker::TRACKING_STATUS::NOT_INIT;
          track_status = false;
        }
        else
        {
          std::cout<<"Tracker: System lost!\n";
          trackingStatus = Abstract_Tracker::TRACKING_STATUS::LOST;
          track_status = false;
        }*/

        track_status = false;
      }

    }

    // Set current as previous frame
    mPrevFrame.swap(mCurrentFrame);
    mCurrentFrame.reset();

    // Return if tracking is ok
    return track_status;
  }


  /// INITIALIZATION

  void resetSystemInitialization()
  {
    std::cout<<"Tracker: Reset system initialization process!\n";
    if (init_ref_frame)
    {
      init_ref_frame.reset();
    }
  }

  void setReferenceSystemInitialization(std::shared_ptr<Frame> & frame)
  {
    // Reset all initialization settings
    resetSystemInitialization();
    // Set current frame as the new reference frame for initialization
    init_ref_frame = frame->share_ptr();
    std::cout<<"Tracker: Set new reference initialization frame!\n";
  }

  // Returns true if current frame is used in the initialization process (also if its put as ref image)
  bool trySystemInitialization
  (
    const image::Image<unsigned char> & ima
  )
  {
    if (trackingStatus == Abstract_Tracker::TRACKING_STATUS::NOT_INIT)
    {
      // Clear all data that could be initialized
      // Set new reference initialization frame
      setReferenceSystemInitialization(mCurrentFrame);
      // Set system to INIT
      trackingStatus = Abstract_Tracker::TRACKING_STATUS::INIT;
      return true;

    }
    else if (trackingStatus == Abstract_Tracker::TRACKING_STATUS::INIT)
    {
      // Check that we have actual reference image
      if (!init_ref_frame)
      {
        // We dont have initialization reference frame (SHOULDNT  HAPPEN)
        // Just in case we reset the system and set the current frame as reference frame
        trackingStatus = Abstract_Tracker::TRACKING_STATUS::NOT_INIT;
        return trySystemInitialization(ima);
      }

      // -------------------
      // -- Match features of two images
      // -------------------
      double startTime = omp_get_wtime();

      Hash_Map<size_t,size_t> putative_matches_ref_cur_idx;
      size_t n_put_matches = matching_AllAll_2D_2D(featureExtractor_, init_ref_frame.get(), mCurrentFrame.get(), putative_matches_ref_cur_idx, init_match_desc_ratio_2);

      double stopTime = omp_get_wtime();
      double secsElapsed = stopTime - startTime; // that's all !
      std::cout<<"Tracker: [Matching All-All 2D-2D] ("<<secsElapsed<<" s)\n";
      std::cout<<"Tracker: [Matching All-All 2D-2D] matches: ("<<n_put_matches<<")\n";

      // Check if enough matches with reference image
      if (n_put_matches < init_track_min_matches)
      {
        // Not enough matches with initialization reference frame
        // Set this frame as initialization reference frame and try with next frame
        std::cout<<"Tracker: Not enough matches ("<<n_put_matches<<") - setting this frame as new reference frame\n";
        setReferenceSystemInitialization(mCurrentFrame);
        trackingStatus = Abstract_Tracker::TRACKING_STATUS::INIT;
        return true;
      }

      // -------------------
      // -- Estimate relative pose between images
      // -------------------
      bool b_motion_est_ok = false;
      // Recovered motion
      Mat3 c_R_r;
      Vec3 c_t_r; // Center of origin of w (reference cam) in c (current cam)
      double c_s_r;
      double AC_reprojection_thresh;
      std::vector<size_t> inliers_M;  // Inliers of an estimated model

      // Intrinsics of the camera
      const IntrinsicBase * cam_intrinsic_ref = init_ref_frame->getCameraIntrinsics();
      const IntrinsicBase * cam_intrinsic_cur = mCurrentFrame->getCameraIntrinsics();

      startTime = omp_get_wtime();

      // If calibrated camera we use Essential matrix otherwise Fundamental
      if (init_ref_frame->getCamCalibrated() && mCurrentFrame->getCamCalibrated())
      {
        b_motion_est_ok = estimateRobustRelativePoseHE(init_ref_frame.get(), mCurrentFrame.get(), putative_matches_ref_cur_idx, init_min_cos_angle_pt, init_track_min_matches, c_R_r, c_t_r, inliers_M, AC_reprojection_thresh);
        // Estimate scale of current camera
        c_s_r = 1.0;
      }
      else
      {
        // TODO: Compute HF
      }

      stopTime = omp_get_wtime();
      secsElapsed = stopTime - startTime; // that's all !
      std::cout<<"Tracker: [Motion] ("<<secsElapsed<<" s)\n";

      if(!b_motion_est_ok)
      {
        std::cout<<"Tracker: Motion estimation failed\n";
        return false;
      }

      // -------------------
      // -- Set initial pose estimations
      // -------------------
      // initial frame
      init_ref_frame->setReferenceFrame_cr_Rts(nullptr,Mat3::Identity(),Vec3::Zero(),1.0);
      init_ref_frame->AC_reprojection_thresh_ = std::max<double>(4,AC_reprojection_thresh);
      // current frame
      mCurrentFrame->setReferenceFrame_cr_Rts(nullptr,c_R_r,c_t_r,c_s_r);
      mCurrentFrame->AC_reprojection_thresh_ = std::max<double>(4,AC_reprojection_thresh);

      std::cout<<"Tracker: [F_1] AC thresh: "<<init_ref_frame->AC_reprojection_thresh_ <<"\n";
      std::cout<<"Tracker: [F_2] AC thresh: "<<mCurrentFrame->AC_reprojection_thresh_ <<"\n";


      // -------------------
      // -- Triangulate inliers
      // -------------------

      // Estimate projection matrix of both cameras based on estimated poses
      Mat34 P1,P2;
      P1 = cartographer_->getCameraProjectionMatrix(init_ref_frame.get(), nullptr);
      P2 = cartographer_->getCameraProjectionMatrix(mCurrentFrame.get(), nullptr);

      // Info about frames
      const size_t frame_init_id = init_ref_frame->getFrameId();
      const size_t frame_cur_id = mCurrentFrame->getFrameId();
      bool cam_init_calibrated = init_ref_frame->getCamCalibrated();
      bool cam_cur_calibrated = mCurrentFrame->getCamCalibrated();

      // Iterate through inliers and triangulate
      std::vector<std::unique_ptr<MapLandmark> > vec_new_triangulated_pts(inliers_M.size());
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for ( size_t k_i = 0; k_i < inliers_M.size(); ++k_i)
      {
        std::unique_ptr<MapLandmark> & ml = vec_new_triangulated_pts[k_i];
        // Initialize new map landmark
        ml.reset(new MapLandmark());
        std::cout<<"NEWT: "<<ml.get()<<"\n";
        // Advance to the correct inlier
        Hash_Map<size_t,size_t>::iterator m_iter = putative_matches_ref_cur_idx.begin();
        std::advance(m_iter,inliers_M[k_i]);

        // Position of detected features
        //   If we calibrated camera we just use undistorted
        //   If not calibrated we remove the distortion as best as we know it (current estimation)
        const size_t & feat_id_ref = m_iter->first;
        const size_t & feat_id_cur = m_iter->second;

        const Vec2
          & x1_ = cam_init_calibrated ? init_ref_frame->getFeaturePosition(feat_id_ref) : cam_intrinsic_ref->remove_disto(init_ref_frame->getFeaturePosition(feat_id_ref)),
          & x2_ = cam_cur_calibrated ? mCurrentFrame->getFeaturePosition(feat_id_cur) : cam_intrinsic_cur->remove_disto(mCurrentFrame->getFeaturePosition(feat_id_cur));

        // Triangulate results
        TriangulateDLT(P1, x1_, P2, x2_, &(ml->X_));

        // Add observations
        LandmarkObservations & obs = ml->obs_;
        obs[frame_init_id] = MapObservation(feat_id_ref,init_ref_frame.get());
        obs[frame_cur_id] = MapObservation(feat_id_cur,mCurrentFrame.get());
      }

      // -------------------
      // -- Optimize initial poses
      // -------------------
      // Set frames that are used in the system
      std::vector<Frame*> local_frames(2);
      local_frames[0] = init_ref_frame.get();
      local_frames[1] = mCurrentFrame.get();

      // HACK -> if we want the first frame to be static we temporarly set it as active
      init_ref_frame->setActive();

      VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
      options.linear_solver_type_ = ceres::DENSE_SCHUR;
      VSSLAM_Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
      if (!bundle_adjustment_obj.OptimizePose(&local_frames,mCurrentFrame.get(), nullptr, &vec_new_triangulated_pts))
      {
        return false;
      }

      // HACK -> if we want the first frame to be static we temporarly set it as active
      init_ref_frame->setInactive();

      // -------------------
      // -- Determine inliers
      // -------------------
      const Mat3 & cam_ref_Rcw = init_ref_frame->getRotationMatrix_rc();
      const Mat3 & cam_cur_Rcw = mCurrentFrame->getRotationMatrix_rc();

      // Go through all the points and check if they are inliers
      for (std::vector<std::unique_ptr<MapLandmark> >::iterator v_it = vec_new_triangulated_pts.begin();
          v_it != vec_new_triangulated_pts.end(); )
      {
        const Vec3 & pt3D = (*v_it)->X_;
        LandmarkObservations & obs = (*v_it)->obs_;

        // Check points for angle and reprojection error
        MapObservation & map_obs_ref = obs.at(frame_init_id);
        MapObservation & map_obs_cur = obs.at(frame_cur_id);

        // Represent points in camera reference frame
        Vec3 pt_3D_ref, pt_3D_cur;
        getRelativePointPosition(pt3D,nullptr,pt_3D_ref,init_ref_frame.get());
        getRelativePointPosition(pt3D,nullptr,pt_3D_cur,mCurrentFrame.get());

        // Check that the point is infront of the cameras
        if (pt_3D_ref(2) <= 0 || pt_3D_cur(2) <= 0)
        {
          v_it = vec_new_triangulated_pts.erase(v_it);
          continue;
        }
        // Get location of features in the image
        const Vec2 & pt_ref = map_obs_ref.frame_ptr->getFeaturePosition(map_obs_ref.feat_id);
        const Vec2 & pt_cur = map_obs_cur.frame_ptr->getFeaturePosition(map_obs_cur.feat_id);

        // Compute rays between points and focal center
        const Vec3 ray_ref = Vec3(cam_ref_Rcw * cam_intrinsic_ref->operator ()(pt_ref)).normalized();
        const Vec3 ray_cur = Vec3(cam_cur_Rcw * cam_intrinsic_cur->operator ()(pt_cur)).normalized();

        // Compute angle between both rays
        const double mag = ray_ref.norm() * ray_cur.norm();
        const double dotAngle = ray_ref.dot( ray_cur );
        const double cosParallax = clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 );

        // If there is not enough parallax between rays we cant robustly triangulate
        if (cosParallax > init_min_cos_angle_pt)
        {
          v_it = vec_new_triangulated_pts.erase(v_it);
          continue;
        }

        // Check reprojection error in reference image
        const Vec2 pt_proj_ref = cam_intrinsic_ref->cam2ima(cam_intrinsic_ref->have_disto()?cam_intrinsic_ref->add_disto(pt_3D_ref.hnormalized()):pt_3D_ref.hnormalized());
        if ((pt_ref - pt_proj_ref).squaredNorm() > init_ref_frame->AC_reprojection_thresh_ )
        {
          v_it = vec_new_triangulated_pts.erase(v_it);
          continue;
        }

        // Check reprojection error in current image
        const Vec2 pt_proj_cur = cam_intrinsic_cur->cam2ima(cam_intrinsic_cur->have_disto()?cam_intrinsic_cur->add_disto(pt_3D_cur.hnormalized()):pt_3D_cur.hnormalized());
        if ((pt_cur - pt_proj_cur).squaredNorm() > mCurrentFrame->AC_reprojection_thresh_ )
        {
          v_it = vec_new_triangulated_pts.erase(v_it);
          continue;
        }

        // Set the best descriptor of the landmark
        cartographer_->updateBestMapPointDescriptor((*v_it).get());

        v_it++;
      }
      std::cout<<"Tracker: Number of init map points: "<<vec_new_triangulated_pts.size()<<"\n";

      // -------------------
      // -- Try to initialize map
      // -------------------
      if (vec_new_triangulated_pts.size() < init_track_min_matches)
      {
        // Not enough matches with initialization reference frame
        // Set this frame as initialization reference frame and try with next frame
        std::cout<<"Tracker: Not enough triangulated points ("<<n_put_matches<<")\n"
            <<" - Setting this frame as new reference frame\n";
        setReferenceSystemInitialization(mCurrentFrame);
        trackingStatus = Abstract_Tracker::TRACKING_STATUS::INIT;
        return true;
      }

      // Add initialization steps for map initialization
      cartographer_->initializationAddStep(init_ref_frame, nullptr);
      cartographer_->initializationAddStep(mCurrentFrame, &vec_new_triangulated_pts);

      // Update motion model
      motionModel.updateMotionModel(init_ref_frame.get(),mCurrentFrame.get());

      // Set second camera as reference keyframe
      mLastRefFrame = mCurrentFrame->share_ptr();

      // Clear initialization data
      resetSystemInitialization();

      // Set tracking to OK
      trackingStatus = Abstract_Tracker::TRACKING_STATUS::OK;
      return true;
    }
  }


  bool trackWithMotionModel()
  {
    Hash_Map<MapLandmark *,size_t> putative_matches_3D_ptr_cur_idx;

    // -------------------
    // -- Predict location of current frame (using MM)
    // -------------------
    Mat4 T_predict = motionModel.predictLocation(mPrevFrame.get(),mCurrentFrame.get());
    mCurrentFrame->setPose_cr_T(T_predict);

    // -------------------
    // -- Match by projecting triangulated points from prev frame to current frame
    // -------------------
    matching_Projection_3D_2D(featureExtractor_, mPrevFrame->map_points_, mCurrentFrame.get(), putative_matches_3D_ptr_cur_idx, track_mm_win_size, track_match_desc_ratio_2);

    std::cout<<"Tracker: Matched by projection: "<<putative_matches_3D_ptr_cur_idx.size()<<"\n";


    // If not enough matches we try again with wider search window
    if (putative_matches_3D_ptr_cur_idx.size() < track_min_matches)
    {
      std::cout<<"Tracker: Matched by projection try 1 failed! Matches: "<<putative_matches_3D_ptr_cur_idx.size()<<"\n";

      // Double search window for searching by projection
      matching_Projection_3D_2D(featureExtractor_, mPrevFrame->map_points_, mCurrentFrame.get(), putative_matches_3D_ptr_cur_idx, track_mm_win_size*2, track_match_desc_ratio_2);
      if (putative_matches_3D_ptr_cur_idx.size() < track_min_matches)
      {
        std::cout<<"Tracker: Matched by projection try 2 failed! Matches: "<<putative_matches_3D_ptr_cur_idx.size()<<"\n";
        return false;
      }
    }

    std::cout<<"Tracker: Matched by projection OK: "<<putative_matches_3D_ptr_cur_idx.size()<<"\n";

    // -------------------
    // -- Do tiny BA to refine the pose based on the matches
    // -------------------
    VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
    VSSLAM_Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    if (!bundle_adjustment_obj.OptimizePose(nullptr, mCurrentFrame.get(), &putative_matches_3D_ptr_cur_idx, nullptr))
    {
      return false;
    }

    // TODO: Check how to set reprojection threshold if we just use BA
    mCurrentFrame->AC_reprojection_thresh_ = default_reproj_thresh_2;

    // Associate the features in frame with the 3D landmarks (if inliers)
    checkReprojectionAndMarkInliersInFrame(mCurrentFrame.get(),putative_matches_3D_ptr_cur_idx);

    size_t n_matches = mCurrentFrame->getNumberMapPoints();

    std::cout<<"Tracker: Matched after BA: "<<n_matches<<"\n";

    if (n_matches < track_min_matches)
      return false;
    return true;
  }

  void checkReprojectionAndMarkInliersInFrame(Frame * frame, Hash_Map<MapLandmark *,size_t> & putative_matches_3D_ptr_cur_idx)
  {
    // -------------------
    // -- Determine which of the points are actually inliers
    // -------------------
    const size_t & frame_cur_id = frame->getFrameId();
    IntrinsicBase * cam_cur_intrinsic = frame->getCameraIntrinsics();
    const double & frame_reproj_thresh = frame->AC_reprojection_thresh_;

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t pm_i = 0; pm_i < putative_matches_3D_ptr_cur_idx.size(); ++pm_i)
    {
      Hash_Map<MapLandmark *,size_t>::iterator iter_p_match = putative_matches_3D_ptr_cur_idx.begin();
      std::advance(iter_p_match,pm_i);
      // Get Landmark data
      MapLandmark * map_point = iter_p_match->first;
      const size_t feat_id_cur = iter_p_match->second;

      // Project point to frame coordinate system
      Vec3 pt_3D_cur;
      getRelativePointPosition(map_point->X_,map_point->ref_frame_,pt_3D_cur,frame);

      // Check that the point is infront of the camera
      if (pt_3D_cur(2) <= 0)
      {
        continue;
      }

      // Get the observation
      const Vec2 & obs_cur = frame->getFeaturePosition(feat_id_cur);

      // Compute residual error in current frame
      // We add distortion to 3D points if we have it
      const Vec2 pt_cur = cam_cur_intrinsic->cam2ima(cam_cur_intrinsic->have_disto()?cam_cur_intrinsic->add_disto(pt_3D_cur.hnormalized()):pt_3D_cur.hnormalized());

      if ((obs_cur - pt_cur).squaredNorm() < frame_reproj_thresh)
      {
        frame->map_points_[feat_id_cur] = map_point;
        map_point->localMapFrameId_ = frame_cur_id;
      }
    }
  }

  void markInliersInFrame(Frame * frame, Hash_Map<MapLandmark *,size_t> & putative_matches_3D_ptr_cur_idx)
  {
    const size_t frame_cur_id = frame->getFrameId();
    // Mark all matches as valid points (they are inside the reprojection error tested in matching_Projection_3D_2D)
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t k_i = 0; k_i < putative_matches_3D_ptr_cur_idx.size(); ++k_i)
    {
      Hash_Map<MapLandmark *, size_t>::iterator it_p_match = putative_matches_3D_ptr_cur_idx.begin();
      std::advance(it_p_match,k_i);

      frame->map_points_[it_p_match->second] = it_p_match->first;
      it_p_match->first->localMapFrameId_ = frame_cur_id;
    }
  }
  bool trackWithReferenceFrame()
  {
    // -------------------
    // -- Set location of current frame (as the one of last frame)
    // -------------------
    mCurrentFrame->setPose_cr_T(mPrevFrame->getTransformationMatrix_cr());

    // -------------------
    // -- Try to match all-all with the last reference frame
    // -------------------
    Hash_Map<MapLandmark* ,size_t> putative_matches_3D_ptr_cur_idx;
    matching_AllAll_3D_2D(featureExtractor_, mLastRefFrame->map_points_, mCurrentFrame.get(), putative_matches_3D_ptr_cur_idx, track_match_desc_ratio_2);

    std::cout<<"Tracker: Matched All-All frames: "<<mLastRefFrame->getFrameId()<<" and "<<mCurrentFrame->getFrameId()<<" :: "<<putative_matches_3D_ptr_cur_idx.size()<<"\n";

    // Check if we got enough matches
    if (putative_matches_3D_ptr_cur_idx.size() < track_min_matches)
    {
      std::cout<<"Tracker: Matched All-All failed! Matches: "<<putative_matches_3D_ptr_cur_idx.size()<<"\n";
      return false;
    }

    // -------------------
    // -- Do tiny BA to refine the pose based on the matches
    // -------------------
    VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
    VSSLAM_Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    if (!bundle_adjustment_obj.OptimizePose(nullptr, mCurrentFrame.get(),&putative_matches_3D_ptr_cur_idx, nullptr))
    {
      return false;
    }

    // TODO: Check how to set reprojection threshold if we just use BA
    mCurrentFrame->AC_reprojection_thresh_ = default_reproj_thresh_2;


    /*
    // -------------------
    // -- Estimate the pose using AC EPnP
    // -------------------
    Mat3 R;
    Vec3 t;
    double s;
    std::vector<size_t> vec_inliers;
    double AC_threshold;

    // TODO: How to estimate scale!!
    bool bSuccessMotion = estimateRobustPose(mCurrentFrame.get(),putative_matches_3D_ptr_cur_idx, track_min_matches,R,t,vec_inliers,AC_threshold);
    s = 1.0f;

    if (!bSuccessMotion)
    {
      std::cout<<"Tracker: Robust pose estimation (AC EPnP) with AllAll failed!\n";
      return false;
    }
    // -------------------
    // -- Update estimated pose
    // -------------------
    mCurrentFrame->setReferenceFrame_cr_Rts(nullptr,R,t,s);
    mCurrentFrame->AC_reprojection_thresh_ = std::max<double>(4,AC_threshold);
    std::cout<<"Tracker: [F_c] AC thresh: "<<mCurrentFrame->AC_reprojection_thresh_ <<"\n";


        // Mark inliers in the set
        const size_t & frame_cur_id = mCurrentFrame->getFrameId();
        #ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        //for (size_t k_i = 0; k_i < putative_matches_3D_ptr_cur_idx.size(); ++k_i)
        for (size_t k_i = 0; k_i < vec_inliers.size(); ++k_i)
        {
          Hash_Map<MapLandmark *, size_t>::iterator it_p_match = putative_matches_3D_ptr_cur_idx.begin();
          //std::advance(it_p_match,k_i);
          std::advance(it_p_match,vec_inliers[k_i]);
          MapLandmark * map_point = it_p_match->first;

          mCurrentFrame->map_points_[it_p_match->second] = map_point;
          map_point->localMapFrameId_ = frame_cur_id;
        }
    */

    // Associate the features in frame with the 3D landmarks (if inliers)
    checkReprojectionAndMarkInliersInFrame(mCurrentFrame.get(),putative_matches_3D_ptr_cur_idx);

    size_t n_matches = mCurrentFrame->getNumberMapPoints();

    std::cout<<"Tracker: Matched after referenceKF: "<<n_matches<<"\n";

    if (n_matches < track_min_matches)
      return false;
    return true;
  }

  void trackLocalMap()
  {
    std::cout<<"Tracker: Track local map\n";

    // -------------------
    // -- Identify local map (frames & points)
    // -------------------
    std::vector<Frame *> local_map_frames;
    std::vector<MapLandmark *> local_map_points;

    // Identify N best frames that see points seen from current camera
    mCurrentFrame->getFrameVisibilityConnections(local_map_frames,track_local_map_size);

    // Get all unmatched map points from neighbor frames
    cartographer_->getLocalMapPoints(mCurrentFrame.get(),local_map_frames,local_map_points);

    std::cout<<"Tracker: Local frames: "<<local_map_frames.size()<<" Local Points: "<<local_map_points.size()<<"\n";

    // -------------------
    // -- Try to match local map points with the features in the current frame
    // -------------------
    Hash_Map<MapLandmark *,size_t> putative_matches_3D_ptr_cur_idx;
    matching_Projection_3D_2D(featureExtractor_, local_map_points, mCurrentFrame.get(), putative_matches_3D_ptr_cur_idx, mCurrentFrame->AC_reprojection_thresh_*2, track_match_desc_ratio_2);

    std::cout<<"Tracker: Matches with local map: "<<putative_matches_3D_ptr_cur_idx.size()<<"\n";

    // Mark all points as inliers
    markInliersInFrame(mCurrentFrame.get(),putative_matches_3D_ptr_cur_idx);

    size_t n_matches = mCurrentFrame->getNumberMapPoints();

    std::cout<<"Tracker: Total matches with map: "<<n_matches<<"\n";
  }


  void triangulateNewPoints(std::vector<Frame *> & local_map_frames, std::vector<std::unique_ptr<MapLandmark> > & vec_new_pts_3D_obs)
  {
    // -------------------
    // -- Identify local map (frames)
    // -------------------
    // Identify N best frames that see points seen from current camera
    mCurrentFrame->getFrameVisibilityConnections(local_map_frames,triangule_local_map_size);

    // -------------------
    // -- Loop through each neighbor frame and check if we can match any of the (unmatched) points in current frame to
    // -- (unmatched) features in the neighbor frames
    // -------------------

    std::vector<Hash_Map<size_t,size_t> > neighbor_matches(local_map_frames.size());

    double startTime = omp_get_wtime();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t f_i = 0; f_i < local_map_frames.size(); ++f_i)
    {
      Frame * frame_neighbor = local_map_frames[f_i];
      // Compute fundamental matrix between F_current and F_neigh
      Mat3 F_n_c; // p_neigh' F_cn p_cur
      computeFundamentalMatrix(mCurrentFrame.get(),frame_neighbor,F_n_c);
      matching_EpipolarLine_2D_2D(featureExtractor_, mCurrentFrame.get(), frame_neighbor, F_n_c, neighbor_matches[f_i],track_epipolar_desc_ratio,featureExtractor_->max_dist_desc_d2);
      std::cout<<"Tracker: EPI match: "<<mCurrentFrame->getFrameId()<<" :: "<<frame_neighbor->getFrameId()<<" :: "<<neighbor_matches[f_i].size()<<"\n";
    }

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Tracker: Neighbour TOTAL ("<<secsElapsed<<")\n";

    // -------------------
    // -- Concatenate observations of the same point over all neighbor frames
    // -------------------
    Hash_Map<size_t,LandmarkObservations> possible_matches_nviews;

    startTime = omp_get_wtime();
    for (size_t f_i = 0; f_i < local_map_frames.size(); ++f_i)
    {
      Frame * frame_i = local_map_frames[f_i];

      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for (size_t m_i=0; m_i < neighbor_matches[f_i].size(); ++m_i)
      {
        Hash_Map<size_t,size_t>::iterator match_i = neighbor_matches[f_i].begin();
        std::advance(match_i,m_i);
        const size_t feat_cur_id = match_i->first;

        if (possible_matches_nviews.find(feat_cur_id)!=possible_matches_nviews.end())
          continue;

        LandmarkObservations obs;
        obs[frame_i->getFrameId()] = MapObservation(match_i->second, frame_i);

        for (size_t f_j = f_i+1; f_j < local_map_frames.size(); ++f_j)
        {
          Hash_Map<size_t,size_t>::iterator it_f = neighbor_matches[f_j].find(feat_cur_id);
          if (it_f != neighbor_matches[f_j].end())
          {
            obs[local_map_frames[f_j]->getFrameId()] = MapObservation(it_f->second, local_map_frames[f_j]);
          }
        }
        if(obs.size()>0)
        #pragma omp critical
        {
          possible_matches_nviews[feat_cur_id] = obs;
        }
      }
    }
    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Tracker: Neighbour filter ("<<secsElapsed<<")\n";

    // -------------------
    // -- Triangulate new points
    // --  Go through measurements and triangulate the first match that has enough angle
    // -------------------
    vec_new_pts_3D_obs.reserve(possible_matches_nviews.size());

    // Estimate projection matrix of both cameras based on estimated poses
    Mat34 Pc,Pi;
    // -------------------
    // -- Loop through matches and decide which we should triangulate
    // -------------------
    Mat3 cam_cur_Rcw = mCurrentFrame->getRotationMatrix_rc();

    const Pinhole_Intrinsic * cam_intrinsic_cur = dynamic_cast<const Pinhole_Intrinsic*>(mCurrentFrame->getCameraIntrinsics());
    const Mat3 K_cur_inv = cam_intrinsic_cur->Kinv();

    // Determine projection matrix of current frame
    //TODO: How to use if we have relative cameras
    Pc = mCurrentFrame->getProjectionMatrix();

    Frame * frame_i;
    Pinhole_Intrinsic * cam_intrinsic_i;
    Mat3 K_i_inv;
    Mat3 cam_i_Rcw;

    size_t t_i = 0;
    for(auto p_matches : possible_matches_nviews)
    {
      const Vec2 & x_cur_ = mCurrentFrame->getCamCalibrated() ? mCurrentFrame->getFeaturePosition(p_matches.first) : cam_intrinsic_cur->remove_disto(mCurrentFrame->getFeaturePosition(p_matches.first));
      // All other views that see the same point
      LandmarkObservations & obs = p_matches.second;

      // Find first pair of images with sufficient angle to triangulate -> add the rest
      for (auto & m_o : obs)
      {
        if (frame_i != m_o.second.frame_ptr)
        {
          frame_i = m_o.second.frame_ptr;
          cam_i_Rcw = frame_i->getRotationMatrix_rc();
          if (cam_intrinsic_i != frame_i->getCameraIntrinsics())
          {
            cam_intrinsic_i = dynamic_cast<Pinhole_Intrinsic*>(frame_i->getCameraIntrinsics());
            K_i_inv = cam_intrinsic_i->Kinv();
          }
          // Get projection matrix of second camera
          //TODO: How to use if we have relative cameras
          Pi = frame_i->getProjectionMatrix();
        }

        // Get measurement in the i-th view
        const Vec2 & x_i_ = frame_i->getCamCalibrated() ? frame_i->getFeaturePosition(m_o.second.feat_id) : cam_intrinsic_i->remove_disto(frame_i->getFeaturePosition(m_o.second.feat_id));

        // Compute ray angle between points
        const Vec3 ray1 = Vec3(cam_cur_Rcw * Vec3( K_cur_inv * Vec3( x_cur_( 0 ), x_cur_( 1 ), 1.0 ) )).normalized();
        const Vec3 ray2 = Vec3(cam_i_Rcw * Vec3( K_i_inv * Vec3( x_i_( 0 ), x_i_( 1 ), 1.0 ) )).normalized();
        const double mag = ray1.norm() * ray2.norm();
        const double dotAngle = ray1.dot( ray2 );
        const double cosParallax = clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 );
        // if too small we dont consider it
        if (cosParallax > init_min_cos_angle_pt)
          continue;

        // Triangulate results
        Vec3 pt_3D;
        TriangulateDLT(Pc, x_cur_, Pi, x_i_, &pt_3D);

        Vec3 pt_3D_cur, pt_3D_i;
        getRelativePointPosition(pt_3D,nullptr,pt_3D_cur,mCurrentFrame.get());
        getRelativePointPosition(pt_3D,nullptr,pt_3D_i,frame_i);

        if (pt_3D_cur(2) <= 0 || pt_3D_i(2) <= 0)
          continue;

        std::unique_ptr<MapLandmark> ml = std::unique_ptr<MapLandmark>(new MapLandmark());
        ml->X_ = pt_3D;
        ml->obs_ = obs;
        ml->localMapFrameId_ = mCurrentFrame->getFrameId();
        // Add current observation
        ml->obs_[mCurrentFrame->getFrameId()] = MapObservation(p_matches.first,mCurrentFrame.get());
        vec_new_pts_3D_obs.push_back(std::move(ml));
        break;
      }

    }
    vec_new_pts_3D_obs.shrink_to_fit();
    std::cout<<"Tracker: Actual new points: "<<vec_new_pts_3D_obs.size()<<"\n";
  }




  void removeOutliersFromFrames(std::vector<Frame *> & local_map_frames, std::vector<std::unique_ptr<MapLandmark> > & vec_new_pts_3D_putative)
  {
    // Go through local frames and check the ones which are not fixed
    for (Frame * & frame_i : local_map_frames)
    {
      // If frame is active the pose is fixed (info in the local ba cant change)
      if (frame_i->isActive())
        continue;

      // -------------------
      // -- Determine which of the points are outliers and remove them
      // -------------------
      const size_t & frame_i_id = frame_i->getFrameId();
      IntrinsicBase * cam_i_intrinsic = frame_i->getCameraIntrinsics();
      const double & frame_i_reproj_error = frame_i->AC_reprojection_thresh_;

      std::vector<MapLandmark*> & frame_i_3D_pts = frame_i->map_points_;

      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for (size_t pm_i = 0; pm_i < frame_i_3D_pts.size(); ++pm_i)
      {
        // Get Landmark data
        MapLandmark * & map_point = frame_i_3D_pts[pm_i];
        if (!map_point)
          continue;

        Vec3 pt_3D_cur;
        getRelativePointPosition(map_point->X_,map_point->ref_frame_,pt_3D_cur,frame_i);

        // Check that the point is infront of the camera
        if (pt_3D_cur(2) < 0)
        {
          // Mark as not matched
          frame_i_3D_pts[pm_i] = nullptr;
          continue;
        }

        // Get the observation
        const Vec2 & obs_cur = frame_i->getFeaturePosition(pm_i);

        // Compute residual error in current frame
        // We add distortion to 3D points if we have it
        const Vec2 pt_cur = cam_i_intrinsic->cam2ima(cam_i_intrinsic->have_disto()?cam_i_intrinsic->add_disto(pt_3D_cur.hnormalized()):pt_3D_cur.hnormalized());

        // Keep if in reprojection threshold
        if ((obs_cur - pt_cur).squaredNorm() < frame_i_reproj_error)
        {
          continue;
        }

        // Mark as not matched
        frame_i_3D_pts[pm_i] = nullptr;
      }
    }

    // Check newly triangulated points if they have at least two inliers
    std::vector<size_t> vec_new_pts_outliers;
    //#ifdef OPENMVG_USE_OPENMP
    //#pragma omp parallel for schedule(dynamic)
    //#endif
    for(size_t t_i = 0; t_i < vec_new_pts_3D_putative.size(); ++t_i)
    {
      MapLandmark * ml = vec_new_pts_3D_putative[t_i].get();

      LandmarkObservations & obs = ml->obs_;

      Vec3 pt_3D_cam;
      for(LandmarkObservations::iterator iter_mo = obs.begin(); iter_mo != obs.end();)
      {
        MapObservation & m_o =  iter_mo->second;
        Frame * & frame_i = m_o.frame_ptr;
        // Put 3D point into coordinate system of frame
        getRelativePointPosition(ml->X_,ml->ref_frame_,pt_3D_cam,frame_i);

        if (frame_i->getSquaredReprojectionError(pt_3D_cam,m_o.feat_id) > frame_i->AC_reprojection_thresh_)
        {
          // Remove measurement from 3D point
          iter_mo = obs.erase(iter_mo);
        }
        else
        {
          ++iter_mo;
        }
      }

      if (obs.size()<2)
      {
        #pragma omp critical
        {
          vec_new_pts_outliers.push_back(t_i);
        }
      }
    }

    // Remove any triangulated landmarks that dont have enough measurements
    std::vector<std::unique_ptr<MapLandmark> >::iterator it_out;
    for (size_t o_i = 0; o_i < vec_new_pts_outliers.size(); ++o_i)
    {
      it_out = vec_new_pts_3D_putative.begin();
      std::advance(it_out,vec_new_pts_outliers[o_i]-o_i);
      vec_new_pts_3D_putative.erase(it_out);
    }
  }



  bool detect
  (
    const image::Image<unsigned char> & ima,
    Frame * frame,
    const size_t & min_count,
    const size_t & max_count
  )
  {
    //double startTime = omp_get_wtime();

    // Detect feature points
    size_t n_feats_detected = featureExtractor_->detect(ima,frame,min_count,max_count);

    //double stopTime = omp_get_wtime();
    //double secsElapsed = stopTime - startTime; // that's all !
    //std::cout<<"Detect time:"<<secsElapsed<<"\n";

    if (!(n_feats_detected > 0))
    {
      return false;
    }

    //startTime = omp_get_wtime();

    // Describe detected features
    featureExtractor_->describe(ima,frame);

    //stopTime = omp_get_wtime();
    //secsElapsed = stopTime - startTime; // that's all !
    //std::cout<<"Describe time:"<<secsElapsed<<"\n";

    // Undistort points

    frame->updateFeaturesData();

    //startTime = omp_get_wtime();
    //stopTime = omp_get_wtime();
    //secsElapsed = stopTime - startTime; // that's all !
    //std::cout<<"Undistort time:"<<secsElapsed<<"\n";

    if (n_feats_detected > 0)
      return true;
    return false;
  }





};

} // namespace VO
} // namespace openMVG
