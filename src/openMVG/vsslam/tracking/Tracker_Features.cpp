// Copyright (c) 2017 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <openMVG/vsslam/tracking/Tracker_Features.hpp>

namespace openMVG  {
namespace VSSLAM  {

  /// Try to track current point set in the provided image
  /// return false when tracking failed (=> to send frame to relocalization)
  bool Tracker_Features::track
  (
    const image::Image<unsigned char> & ima,
    std::shared_ptr<Frame> current_frame
  )
  {
    bool track_status = false;
    // Set current frame
    mCurrentFrame = current_frame->share_ptr();

    vec_times[0] = omp_get_wtime();
    std::cout<<"Tracker: Start frame: "<<current_frame->getFrameId()<<"\n";
    // -------------------------
    // -- Detect features
    // -------------------------
    if (!detect(ima,mCurrentFrame.get(),max_tracked_points,0))
    {
      // No features detected on the frame
      std::cout<<"Tracker: Error - no features detected!\n";
      return false;
    }
    const size_t & n_feats_detected = mCurrentFrame->getNumberOfFeatures();
    std::cout<<"Tracker: Features detected: "<<n_feats_detected<<"\n";

    std::cout<<"Tracker: [Detection] ("<<omp_get_wtime() - vec_times[0]<<" s)\n";
    vec_times[1] = omp_get_wtime();

    // -------------------------
    // -- Tracking
    // -------------------------
    // If tracking is not initialized
    if (trackingStatus == TRACKING_STATUS::NOT_INIT)
    {
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
      bool b_track_OK = false;  // Flag to determine if we are able to track one way or another

      // ----------------
      // -- Motion Model tracking
      // ----------------
      // Use previous tracking to predict where the pose will be and try matching with previously triangulated map points
      if (motionModel.isValid())
      {
        // Match current frame with lastFrame
        b_track_OK = trackWithMotionModel();
      }
      std::cout<<"Tracker: [MM] ("<<omp_get_wtime() - vec_times[1]<<" s)\n";
      vec_times[1] = omp_get_wtime();

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

      std::cout<<"Tracker: [RF] ("<<omp_get_wtime() - vec_times[1]<<" s)\n";
      vec_times[1] = omp_get_wtime();

      // If sucessful triangulate new points
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
        triangulateNewLandmarks(vec_new_pts_3D);


        std::cout<<"Tracker: [Triangulate new pts] ("<<omp_get_wtime() - vec_times[1]<<" s)\n";
        vec_times[1] = omp_get_wtime();

        // ----------------
        // -- Perform BA on all the matches existing points and new trinagulated points
        // ----------------
        if (!cartographer_->optimizeLocal(mCurrentFrame.get(), vec_new_pts_3D))
        {
          return false;
        }

        std::cout<<"Tracker: [Local optimization] ("<<omp_get_wtime() - vec_times[1]<<" s)\n";
        vec_times[1] = omp_get_wtime();

        // ----------------
        // -- Check local map
        // ----------------
        cartographer_->verifyLocalLandmarks(mCurrentFrame.get());

        //removeOutliersInLocalFrames(local_map_frames);
        std::cout<<"Tracker: [Remove outliers - frames] ("<<omp_get_wtime() - vec_times[1]<<" s)\n";

        vec_times[1] = omp_get_wtime();

        removeOutliersInNewTriangulatedPoints(vec_new_pts_3D);
        std::cout<<"Tracker: [Remove outliers - triangulated points] ("<<omp_get_wtime() - vec_times[1]<<" s)\n";
        vec_times[1] = omp_get_wtime();

        std::cout<<"Tracker: New Trianguulated points: "<<vec_new_pts_3D.size()<<"\n";


        size_t n_pts_frame = mCurrentFrame->getNumberMapPoints();
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

            std::cout<<"Tracker: [Init ADD] ("<<omp_get_wtime() - vec_times[1]<<" s)\n";
            vec_times[1] = omp_get_wtime();
          }
          else
          {
            cartographer_->addStep(mCurrentFrame, &vec_new_pts_3D);
          }


          std::cout<<"Tracker: [ADD] ("<<omp_get_wtime() - vec_times[1]<<" s)\n";

          // Export step to Ply
          std::ostringstream os;
          os << std::setw(8) << std::setfill('0') << "Scene_"<<cartographer_->step_id;
          cartographer_->exportSceneToPly(stlplus::create_filespec("/home/klemen/exp_ply", os.str(), ".ply"), sfm::ESfM_Data(sfm::ALL));


          vec_times[1] = omp_get_wtime();
          //mCurrentFrame->computeFrameVisibilityConnections();

          std::cout<<"Tracker: [Compute Visibility] ("<<omp_get_wtime() - vec_times[1]<<" s)\n";
          vec_times[1] = omp_get_wtime();
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
    std::cout<<"Tracker: [End tracker] ("<<omp_get_wtime() - vec_times[1]<<" s)\n";
    vec_times[1] = omp_get_wtime();

    // Return if tracking is ok
    return track_status;
  }


  /// INITIALIZATION
  void Tracker_Features::resetSystemInitialization()
  {
    std::cout<<"Tracker: Reset system initialization process!\n";
    if (init_ref_frame)
    {
      init_ref_frame.reset();
    }
  }

  void Tracker_Features::setReferenceSystemInitialization(std::shared_ptr<Frame> & frame)
  {
    // Reset all initialization settings
    resetSystemInitialization();
    // Set current frame as the new reference frame for initialization
    init_ref_frame = frame->share_ptr();
    std::cout<<"Tracker: Set new reference initialization frame!\n";
  }

  // Returns true if current frame is used in the initialization process (also if its put as ref image)
  bool Tracker_Features::trySystemInitialization
  (
    const image::Image<unsigned char> & ima
  )
  {
    if (trackingStatus == Abstract_Tracker::TRACKING_STATUS::NOT_INIT)
    {
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

      double start_time = omp_get_wtime();
      // -------------------
      // -- Match features of two images
      // -------------------
      matching::IndMatches vec_putative_matches_ref_cur_idx;
      // Match all-all features
      featureMatcher_->matching_AllAll_2D_2D(featureExtractor_, init_ref_frame.get(), mCurrentFrame.get(), vec_putative_matches_ref_cur_idx, init_match_desc_ratio);
      size_t n_put_matches = vec_putative_matches_ref_cur_idx.size();
      std::cout<<"Tracker: [Matching All-All 2D-2D] matches: ("<<n_put_matches<<")\n";

      std::cout<<"Tracker: [Matching All-All 2D-2D] ("<<omp_get_wtime() - start_time<<" s)\n";

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

      // Check if cameras are calibrated
      bool b_cam_init_calibrated = init_ref_frame->getCamCalibrated();
      bool b_cam_cur_calibrated = mCurrentFrame->getCamCalibrated();


      start_time = omp_get_wtime();
      // Recovered motion
      Mat4 c_T_r; // current_T_reference
      double AC_reproj_thresh_2;  // Reprojection error estimated for the model
      std::vector<size_t> inliers_M;  // Inliers of an estimated model

      // If calibrated camera we use Essential matrix otherwise Fundamental
      if (b_cam_init_calibrated && b_cam_cur_calibrated)
      {
        b_motion_est_ok = estimateRobustRelativePosePinholeHE(init_ref_frame.get(), mCurrentFrame.get(), vec_putative_matches_ref_cur_idx, init_min_cos_angle_pt, init_track_min_matches, c_T_r, inliers_M, AC_reproj_thresh_2);
      }
      else
      {
        // TODO: Compute HF
      }

      std::cout<<"Tracker: [Motion] ("<<omp_get_wtime() - start_time<<" s)\n";


      if(!b_motion_est_ok)
      {
        std::cout<<"Tracker: Motion estimation failed\n";
        return false;
      }

      // -------------------
      // -- Set initial pose estimations
      // -------------------
      start_time = omp_get_wtime();
      // initial frame
      init_ref_frame->setReferenceFrame_cr_T(nullptr,Mat4::Identity());
      init_ref_frame->AC_reprojection_thresh_ = std::max<double>(4,AC_reproj_thresh_2);
      // current frame
      mCurrentFrame->setReferenceFrame_cr_T(nullptr,c_T_r);
      mCurrentFrame->AC_reprojection_thresh_ = std::max<double>(4,AC_reproj_thresh_2);

      std::cout<<"Tracker: [F_1] AC thresh: "<<init_ref_frame->AC_reprojection_thresh_ <<"\n";
      std::cout<<"Tracker: [F_2] AC thresh: "<<mCurrentFrame->AC_reprojection_thresh_ <<"\n";

      std::cout<<"Tracker: [Init Pose] ("<<omp_get_wtime() - start_time<<" s)\n";

      // -------------------
      // -- Triangulate inliers
      // -------------------
      start_time = omp_get_wtime();

      // Info about frames
      const IndexT frame_init_id = init_ref_frame->getFrameId();
      const IndexT frame_cur_id = mCurrentFrame->getFrameId();

      // Intrinsics of the camera
      const IntrinsicBase * cam_intrinsic_ref = init_ref_frame->getCameraIntrinsics();
      const IntrinsicBase * cam_intrinsic_cur = mCurrentFrame->getCameraIntrinsics();

      // Estimate projection matrix of both cameras based on estimated poses
      Mat34 P1,P2;
      P1 = cartographer_->getCameraProjectionMatrix(init_ref_frame.get(), nullptr);
      P2 = cartographer_->getCameraProjectionMatrix(mCurrentFrame.get(), nullptr);



      // Iterate through inliers and triangulate
      std::vector<std::unique_ptr<MapLandmark> > vec_new_triangulated_pts(inliers_M.size());
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for ( size_t k_i = 0; k_i < inliers_M.size(); ++k_i)
      {
        // Advance to the correct inlier
        matching::IndMatches::iterator m_iter = vec_putative_matches_ref_cur_idx.begin();
        std::advance(m_iter,inliers_M[k_i]);
        // Inlier feat ids
        const IndexT & feat_id_ref = m_iter->i_;
        const IndexT & feat_id_cur = m_iter->j_;

        // Initialize new map landmark
        std::unique_ptr<MapLandmark> & ml = vec_new_triangulated_pts[k_i];
        ml.reset(new MapLandmark());

        // Position of detected features
        //   If we calibrated camera we just use undistorted
        //   If not calibrated we remove the distortion as best as we know it (current estimation)
        const Vec2
          & pt_1 = b_cam_init_calibrated ? init_ref_frame->getFeaturePosition(feat_id_ref) : cam_intrinsic_ref->remove_disto(init_ref_frame->getFeaturePosition(feat_id_ref)),
          & pt_2 = b_cam_cur_calibrated ? mCurrentFrame->getFeaturePosition(feat_id_cur) : cam_intrinsic_cur->remove_disto(mCurrentFrame->getFeaturePosition(feat_id_cur));

        // Triangulate results
        TriangulateDLT(P1, pt_1, P2, pt_2, &(ml->X_));

        // Add observations
        LandmarkObservations & obs = ml->obs_;
        obs[frame_init_id] = MapObservation(feat_id_ref,init_ref_frame.get());
        obs[frame_cur_id] = MapObservation(feat_id_cur,mCurrentFrame.get());

        // Mark as new initialization pt
        ml->association_type_ = 1;
      }
      std::cout<<"Tracker: [Init Triangulation] ("<<omp_get_wtime() - start_time<<" s)\n";

      // -------------------
      // -- Optimize initial poses
      // -------------------
      // First frame is fixed
      init_ref_frame->setActive();

      start_time = omp_get_wtime();
      if (!cartographer_->optimizeLocal(mCurrentFrame.get(), vec_new_triangulated_pts))
      {
        return false;
      }

      std::cout<<"Tracker: [Optimize pose] ("<<omp_get_wtime() - start_time<<" s)\n";

      // -------------------
      // -- Determine inliers
      // -------------------
      start_time = omp_get_wtime();

      removeOutliersInNewTriangulatedPoints(vec_new_triangulated_pts);

      std::cout<<"Tracker: [Outlier detection] ("<<omp_get_wtime() - start_time<<" s)\n";

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
      start_time = omp_get_wtime();
      cartographer_->initializationAddStep(init_ref_frame, nullptr);
      cartographer_->initializationAddStep(mCurrentFrame, &vec_new_triangulated_pts);

      std::cout<<"Tracker: [Map init 1&2] ("<<omp_get_wtime() - start_time<<" s)\n";

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


  bool Tracker_Features::trackWithMotionModel()
  {
    Hash_Map<MapLandmark *,IndexT> map_putative_matches_3D_ptr_cur_idx;

    // -------------------
    // -- Predict location of current frame (using MM)
    // -------------------
    Mat4 T_predict = motionModel.predictLocation(mPrevFrame.get(),mCurrentFrame.get());
    mCurrentFrame->setPose_cr_T(T_predict);

    // -------------------
    // -- Match by projecting triangulated points from prev frame to current frame
    // -------------------
    double start_time = omp_get_wtime();
    featureMatcher_->matching_Projection_3D_2D(featureExtractor_, mPrevFrame->map_points_, mCurrentFrame.get(), map_putative_matches_3D_ptr_cur_idx, track_mm_win_size, track_match_desc_ratio, featureExtractor_->max_dist_desc_);
    std::cout<<"Tracker: [Match Projection 3D-2D - Test 1] ("<<omp_get_wtime() - start_time<<" s)\n";

    // If not enough matches we try again with wider search window
    if (map_putative_matches_3D_ptr_cur_idx.size() < track_min_matches)
    {
      std::cout<<"Tracker: Matched by projection try 1 failed! Matches: "<<map_putative_matches_3D_ptr_cur_idx.size()<<"\n";

      // Double search window for searching by projection
      start_time = omp_get_wtime();
      featureMatcher_->matching_Projection_3D_2D(featureExtractor_, mPrevFrame->map_points_, mCurrentFrame.get(), map_putative_matches_3D_ptr_cur_idx, track_mm_win_size*4, track_match_desc_ratio, featureExtractor_->max_dist_desc_);

      std::cout<<"Tracker: [Match Projection 3D-2D - Test 2] ("<<omp_get_wtime() - start_time<<" s)\n";

      if (map_putative_matches_3D_ptr_cur_idx.size() < track_min_matches)
      {
        std::cout<<"Tracker: Matched by projection try 2 failed! Matches: "<<map_putative_matches_3D_ptr_cur_idx.size()<<"\n";
        return false;
      }
    }

    std::cout<<"Tracker: Matched by projection OK: "<<map_putative_matches_3D_ptr_cur_idx.size()<<"\n";

    // -------------------
    // -- Do tiny BA to refine the pose based on the matches
    // -------------------

    start_time = omp_get_wtime();
    if (!cartographer_->optimizePose(mCurrentFrame.get(), map_putative_matches_3D_ptr_cur_idx))
    {
      return false;
    }

    std::cout<<"Tracker: [Local optimization] ("<<omp_get_wtime() - start_time<<" s)\n";

    // TODO: Check how to set reprojection threshold if we just use BA
    mCurrentFrame->AC_reprojection_thresh_ = default_reproj_thresh_2;

    // Associate the features in frame with the 3D landmarks (if inliers)

    start_time = omp_get_wtime();
    checkReprojectionAndMarkInliersInFrame(mCurrentFrame.get(),map_putative_matches_3D_ptr_cur_idx);

    std::cout<<"Tracker: [Check inliers] ("<<omp_get_wtime() - start_time<<" s)\n";

    size_t n_matches = mCurrentFrame->getNumberMapPoints();

    std::cout<<"Tracker: Matched after BA: "<<n_matches<<"\n";

    if (n_matches < track_min_matches)
      return false;
    return true;
  }

  void Tracker_Features::checkReprojectionAndMarkInliersInFrame(Frame * frame, Hash_Map<MapLandmark *,IndexT> & map_putative_matches_3D_pts_frame_idx)
  {
    // -------------------
    // -- Determine which of the points are actually inliers
    // -------------------
    const IndexT & frame_id = frame->getFrameId();
    IntrinsicBase * cam_intrinsic = frame->getCameraIntrinsics();
    const double & frame_reproj_thresh = frame->AC_reprojection_thresh_;

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t pm_i = 0; pm_i < map_putative_matches_3D_pts_frame_idx.size(); ++pm_i)
    {
      Hash_Map<MapLandmark *,IndexT>::iterator iter_p_match = map_putative_matches_3D_pts_frame_idx.begin();
      std::advance(iter_p_match,pm_i);
      // Get Landmark data
      MapLandmark * map_point = iter_p_match->first;
      const IndexT feat_id_cur = iter_p_match->second;

      // Project point to frame coordinate system
      Vec3 pt_3D_frame;
      getRelativePointPosition(map_point->X_,map_point->ref_frame_,pt_3D_frame,frame);

      // Check that the point is infront of the camera
      if (pt_3D_frame(2) <= 0)
      {
        continue;
      }

      // Get the observation
      const Vec2 & pt_frame = frame->getFeaturePosition(feat_id_cur);

      // Compute residual error in current frame
      // We add distortion to 3D points if we have it
      const Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(cam_intrinsic->have_disto()?cam_intrinsic->add_disto(pt_3D_frame.hnormalized()):pt_3D_frame.hnormalized());

      if ((pt_frame - pt_3D_frame_projected).squaredNorm() < frame_reproj_thresh)
      {
        frame->setMapPoint(feat_id_cur,map_point);
        // Mark as motion model point or reference keyframe
        map_point->association_type_ = 2;

        //map_point->obs_[frame_id] = MapObservation(feat_id_cur,frame);

      }
    }
  }

  void Tracker_Features::markInliersInFrame(Frame * frame, Hash_Map<MapLandmark *,IndexT> & map_putative_matches_3D_pts_frame_idx)
  {
    const IndexT frame_id = frame->getFrameId();
    // Mark all matches as valid points (they are inside the reprojection error tested in matching_Projection_3D_2D)
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t k_i = 0; k_i < map_putative_matches_3D_pts_frame_idx.size(); ++k_i)
    {
      Hash_Map<MapLandmark *, IndexT>::iterator it_p_match = map_putative_matches_3D_pts_frame_idx.begin();
      std::advance(it_p_match,k_i);

      frame->setMapPoint(it_p_match->second,it_p_match->first);
      //frame->map_points_[it_p_match->second] = it_p_match->first;
      //it_p_match->first->localMapFrameId_ = frame_id;

      // Mark as map tracking point
      it_p_match->first->association_type_ = 3;

      //it_p_match->first->obs_[frame_id] = MapObservation(it_p_match->second,frame);
    }
  }

  bool Tracker_Features::trackWithReferenceFrame()
  {
    // -------------------
    // -- Set location of current frame (as the one of last frame)
    // -------------------
    mCurrentFrame->setPose_cr_T(mPrevFrame->getTransformationMatrix_cr());

    double start_time = omp_get_wtime();
    // -------------------
    // -- Try to match all-all with the last reference frame
    // -------------------
    Hash_Map<MapLandmark* ,IndexT> map_putative_matches_3D_pts_frame_idx;
    featureMatcher_->matching_AllAll_3D_2D(featureExtractor_, mLastRefFrame->map_points_, mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx, track_match_desc_ratio, featureExtractor_->max_dist_desc_);

    std::cout<<"Tracker: Matched All-All frames: "<<mLastRefFrame->getFrameId()<<" (pts: "<<mLastRefFrame->getNumberMapPoints()<<") and "<<mCurrentFrame->getFrameId()<<" :: "<<map_putative_matches_3D_pts_frame_idx.size()<<"\n";

    std::cout<<"Tracker: [Match All-All 3D-2D] ("<<omp_get_wtime() - start_time<<" s)\n";

    // Check if we got enough matches
    if (map_putative_matches_3D_pts_frame_idx.size() < track_min_matches)
    {
      std::cout<<"Tracker: Matched All-All failed! Matches: "<<map_putative_matches_3D_pts_frame_idx.size()<<"\n";
      return false;
    }

    // -------------------
    // -- Do tiny BA to refine the pose based on the matches
    // -------------------
    start_time = omp_get_wtime();
    if (!cartographer_->optimizePose(mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx))
    {
      return false;
    }

    std::cout<<"Tracker: [local Optimize] ("<<omp_get_wtime() - start_time<<" s)\n";

    // TODO: Check how to set reprojection threshold if we just use BA
    mCurrentFrame->AC_reprojection_thresh_ = default_reproj_thresh_2;

    // Associate the features in frame with the 3D landmarks (if inliers)
    start_time = omp_get_wtime();
    checkReprojectionAndMarkInliersInFrame(mCurrentFrame.get(),map_putative_matches_3D_pts_frame_idx);

    std::cout<<"Tracker: [Check Inliers] ("<<omp_get_wtime() - start_time<<" s)\n";


    size_t n_matches = mCurrentFrame->getNumberMapPoints();
    std::cout<<"Tracker: Matched after referenceKF: "<<n_matches<<"\n";

    if (n_matches < track_min_matches)
      return false;
    return true;




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


  }

  void Tracker_Features::trackLocalMap()
  {
    std::cout<<"Tracker: [Track local map]\n";

    double start_time = omp_get_wtime();
    // -------------------
    // -- Identify local map (frames & points)
    // -------------------
    std::vector<Frame *> local_map_frames;
    std::vector<MapLandmark *> local_map_points;

    // Identify N best frames that see points seen from current camera
    mCurrentFrame->getFrameVisibilityConnections(local_map_frames,track_local_map_size);

    std::cout<<"Tracker: [Get visibility] ("<<omp_get_wtime() - start_time<<" s)\n";

    // Get all unmatched map points from neighbor frames
    start_time = omp_get_wtime();
    cartographer_->getLocalMapPoints(mCurrentFrame.get(),local_map_frames,local_map_points);
    std::cout<<"Tracker: Local frames: "<<local_map_frames.size()<<" Local Points: "<<local_map_points.size()<<"\n";

    std::cout<<"Tracker: [Get LP] ("<<omp_get_wtime() - start_time<<" s)\n";

    // -------------------
    // -- Try to match local map points with the features in the current frame
    // -------------------
    Hash_Map<MapLandmark *,IndexT> map_putative_matches_3D_pts_frame_idx;

    start_time = omp_get_wtime();
    featureMatcher_->matching_Projection_3D_2D(featureExtractor_, local_map_points, mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx, track_local_map_win_size, track_match_desc_ratio, featureExtractor_->max_dist_desc_);
    std::cout<<"Tracker: [Match Projection 3D-2D] ("<<omp_get_wtime() - start_time<<" s)\n";

    std::cout<<"Tracker: Matches with local map: "<<map_putative_matches_3D_pts_frame_idx.size()<<"\n";

    // Mark all points as inliers
    start_time = omp_get_wtime();
    markInliersInFrame(mCurrentFrame.get(),map_putative_matches_3D_pts_frame_idx);

    std::cout<<"Tracker: [Mark Inliers] ("<<omp_get_wtime() - start_time<<" s)\n";

    size_t n_matches = mCurrentFrame->getNumberMapPoints();

    std::cout<<"Tracker: Total matches with map: "<<n_matches<<"\n";
  }

  void Tracker_Features::triangulateNewLandmarks(std::vector<std::unique_ptr<MapLandmark> > & vec_new_pts_3D)
  {
    std::vector<Frame *> local_map_frames;

    // -------------------
    // -- Identify local map (frames)
    // -------------------
    // Identify N best frames that see points seen from current camera
    double start_time = omp_get_wtime();
    mCurrentFrame->getFrameVisibilityConnections(local_map_frames,triangule_local_map_size);

    std::cout<<"Tracker: [Get visibility] ("<<omp_get_wtime() - start_time<<" s)\n";

    // -------------------
    // -- Loop through each neighbor frame and check if we can match any of the (unmatched) points in current frame to
    // -- (unmatched) features in the neighbor frames
    // -------------------
    std::vector<Hash_Map<IndexT,IndexT> > vec_neighbor_matches_cur_local(local_map_frames.size());

    start_time = omp_get_wtime();
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t f_i = 0; f_i < local_map_frames.size(); ++f_i)
    {
      Frame * local_frame_i = local_map_frames[f_i];
      // Compute fundamental matrix between F_current and F_neigh
      Mat3 F_l_c; // p_local' F_cn p_cur
      computeFundamentalMatrix(mCurrentFrame.get(),local_frame_i,F_l_c);
      featureMatcher_->matching_EpipolarLine_2D_2D(featureExtractor_, mCurrentFrame.get(), local_frame_i, F_l_c, vec_neighbor_matches_cur_local[f_i],track_epipolar_d2, track_epipolar_desc_ratio);
      std::cout<<"Tracker: EPI match: "<<mCurrentFrame->getFrameId()<<" :: "<<local_frame_i->getFrameId()<<" :: "<<vec_neighbor_matches_cur_local[f_i].size()<<"\n";
    }

    std::cout<<"Tracker: [Match Epipolar 2D 2D] ("<<omp_get_wtime() - start_time<<" s)\n";

    // -------------------
    // -- Concatenate observations of the same point over all neighbor frames
    // -------------------
    start_time = omp_get_wtime();
    Hash_Map<IndexT,LandmarkObservations> map_matches_cur_local_nviews;

    for (size_t local_frame_i = 0; local_frame_i < local_map_frames.size(); ++local_frame_i)
    {
      Frame * local_frame = local_map_frames[local_frame_i];
      const IndexT local_frame_id = local_frame->getFrameId();

      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for (size_t m_i = 0; m_i < vec_neighbor_matches_cur_local[local_frame_i].size(); ++m_i)
      {
        Hash_Map<IndexT,IndexT>::iterator match_i = vec_neighbor_matches_cur_local[local_frame_i].begin();
        std::advance(match_i,m_i);
        const IndexT feat_cur_id = match_i->first;

        if (map_matches_cur_local_nviews.find(feat_cur_id)!=map_matches_cur_local_nviews.end())
          continue;

        LandmarkObservations feat_cur_obs_nviews;
        feat_cur_obs_nviews[local_frame_id] = MapObservation(match_i->second, local_frame);

        for (size_t local_frame_j = local_frame_i+1; local_frame_j < local_map_frames.size(); ++local_frame_j)
        {
          Hash_Map<IndexT,IndexT>::iterator it_matches_cur_frame_j = vec_neighbor_matches_cur_local[local_frame_j].find(feat_cur_id);
          if (it_matches_cur_frame_j != vec_neighbor_matches_cur_local[local_frame_j].end())
          {
            feat_cur_obs_nviews[local_map_frames[local_frame_j]->getFrameId()] = MapObservation(it_matches_cur_frame_j->second, local_map_frames[local_frame_j]);
          }
        }
        if(feat_cur_obs_nviews.size()>0)
        {
          #pragma omp critical
          {
            map_matches_cur_local_nviews[feat_cur_id] = feat_cur_obs_nviews;
          }
        }
      }
    }
    std::cout<<"Tracker: [Neighbour Concatenate] ("<<omp_get_wtime() - start_time<<" s)\n";

    // -------------------
    // -- Triangulate new points
    // --  Go through measurements and triangulate the first match that has enough angle
    // -------------------
    start_time = omp_get_wtime();
    vec_new_pts_3D.reserve(map_matches_cur_local_nviews.size());

    // Estimate projection matrix of both cameras based on estimated poses
    Mat34 P_cur,P_i;
    // -------------------
    // -- Loop through matches and decide which we should triangulate
    // -------------------
    Mat3 frame_cur_Rwc = mCurrentFrame->getRotationMatrix_rc();

    const IntrinsicBase * cam_intrinsic_cur = mCurrentFrame->getCameraIntrinsics();
    const Mat3 K_inv_cur = mCurrentFrame->getK_inv();

    // Determine projection matrix of current frame
    //TODO: How to use if we have relative cameras
    P_cur = mCurrentFrame->getProjectionMatrix();

    Frame * frame_i;
    IntrinsicBase * cam_intrinsic_i;
    Mat3 K_inv_i;
    Mat3 frame_i_Rwc;

    for(auto p_matches : map_matches_cur_local_nviews)
    {
      const Vec2 & pt_cur = mCurrentFrame->getCamCalibrated() ? mCurrentFrame->getFeaturePosition(p_matches.first) : cam_intrinsic_cur->remove_disto(mCurrentFrame->getFeaturePosition(p_matches.first));
      // All other views that see the same point
      LandmarkObservations & pt_cur_obs = p_matches.second;

      // Find first pair of images with sufficient angle to triangulate -> add the rest
      for (LandmarkObservations::iterator it_pt_cur_obs = pt_cur_obs.begin(); it_pt_cur_obs != pt_cur_obs.end();)
      {
        MapObservation & m_o = it_pt_cur_obs->second;
        if (frame_i != m_o.frame_ptr)
        {
          frame_i = m_o.frame_ptr;
          frame_i_Rwc = frame_i->getRotationMatrix_rc();
          if (cam_intrinsic_i != frame_i->getCameraIntrinsics())
          {
            cam_intrinsic_i = frame_i->getCameraIntrinsics();
            K_inv_i = frame_i->getK_inv();
          }
          // Get projection matrix of second camera
          //TODO: How to use if we have relative cameras
          P_i = frame_i->getProjectionMatrix();
        }

        // Get measurement in the i-th view
        const Vec2 & pt_i = frame_i->getCamCalibrated() ? frame_i->getFeaturePosition(m_o.feat_id) : cam_intrinsic_i->remove_disto(frame_i->getFeaturePosition(m_o.feat_id));

        // Compute ray angle between points
        const Vec3 ray_cur = Vec3(frame_cur_Rwc * Vec3( K_inv_cur * Vec3( pt_cur( 0 ), pt_cur( 1 ), 1.0 ) )).normalized();
        const Vec3 ray_i = Vec3(frame_i_Rwc * Vec3( K_inv_i * Vec3( pt_i( 0 ), pt_i( 1 ), 1.0 ) )).normalized();
        const double mag = ray_cur.norm() * ray_i.norm();
        const double dotAngle = ray_cur.dot( ray_i );
        const double cosParallax = clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 );
        // if parallax between two points is too small we dont consider it
        if (cosParallax > init_min_cos_angle_pt)
        {
          it_pt_cur_obs = pt_cur_obs.erase(it_pt_cur_obs);
        }
        else
          ++it_pt_cur_obs;
      }
      // If we dont have at least two valid observations we dont even test the point
      if (pt_cur_obs.empty())
        continue;

      // Triangulate from first two views -> check for the rest if it fits the reprojection error
      MapObservation & m_o = pt_cur_obs.begin()->second;
      if (frame_i != m_o.frame_ptr)
      {
        frame_i = m_o.frame_ptr;
        if (cam_intrinsic_i != frame_i->getCameraIntrinsics())
        {
          cam_intrinsic_i = frame_i->getCameraIntrinsics();
        }
        // Get projection matrix of second camera
        //TODO: How to use if we have relative cameras
        P_i = frame_i->getProjectionMatrix();
      }

      const Vec2 & pt_i = frame_i->getCamCalibrated() ? frame_i->getFeaturePosition(m_o.feat_id) : cam_intrinsic_i->remove_disto(frame_i->getFeaturePosition(m_o.feat_id));

      // Triangulate results
      Vec3 pt_3D;
      TriangulateDLT(P_cur, pt_cur, P_i, pt_i, &pt_3D);

      Vec3 pt_3D_i;
      getRelativePointPosition(pt_3D,nullptr,pt_3D_i,mCurrentFrame.get());
      if (pt_3D_i(2) <= 0)
        continue;

      getRelativePointPosition(pt_3D,nullptr,pt_3D_i,frame_i);
      if (pt_3D_i(2) <= 0)
        continue;

      std::unique_ptr<MapLandmark> ml = std::unique_ptr<MapLandmark>(new MapLandmark());
      ml->X_ = pt_3D;
      ml->obs_ = pt_cur_obs;
      ml->obs_[mCurrentFrame->getFrameId()] = MapObservation(p_matches.first,mCurrentFrame.get());
      ml->last_local_map_frame_id_ = mCurrentFrame->getFrameId();
      // Mark as new triangulated point
      ml->association_type_ = 4;
      // Add current observation
      vec_new_pts_3D.push_back(std::move(ml));
    }
    vec_new_pts_3D.shrink_to_fit();
    std::cout<<"Tracker: Actual new points: "<<vec_new_pts_3D.size()<<"\n";

    std::cout<<"Tracker: [New Pts triangulated] ("<<omp_get_wtime() - start_time<<" s)\n";
  }
/*
  void Tracker_Features::removeOutliersInLocalFrames()
  {
    // Go through local frames and check the ones which are not fixed
    for (auto & it_landmark : tmp_structure)
    {
      Frame * frame_i = it_frame.first;

      // -------------------
      // -- Determine which of the points are outliers and remove them
      // -------------------
      const IndexT & frame_i_id = frame_i->getFrameId();
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

        Vec3 pt_3D_frame_i;
        getRelativePointPosition(map_point->X_,map_point->ref_frame_,pt_3D_frame_i,frame_i);

        // Check that the point is infront of the camera
        if (pt_3D_frame_i(2) <= 0)
        {
          // Mark as not matched
          frame_i_3D_pts[pm_i] = nullptr;
          continue;
        }

        // Get the observation
        const Vec2 & pt_frame_i = frame_i->getFeaturePosition(pm_i);

        // Compute residual error in current frame
        // We add distortion to 3D points if we have it
        const Vec2 pt_3D_frame_i_projected = cam_i_intrinsic->cam2ima(cam_i_intrinsic->have_disto()?cam_i_intrinsic->add_disto(pt_3D_frame_i.hnormalized()):pt_3D_frame_i.hnormalized());

        // Keep if in reprojection threshold
        if ((pt_frame_i - pt_3D_frame_i_projected).squaredNorm() > frame_i_reproj_error)
        {
          // Mark as not matched
          frame_i_3D_pts[pm_i] = nullptr;
        }
      }
    }
  }
*/
  void Tracker_Features::removeOutliersInNewTriangulatedPoints(std::vector<std::unique_ptr<MapLandmark> > & vec_putative_new_pts_3D)
  {
    // Check newly triangulated points if they have at least two inliers
    std::vector<size_t> vec_new_pts_outliers;
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for(size_t t_i = 0; t_i < vec_putative_new_pts_3D.size(); ++t_i)
    {
      MapLandmark * ml = vec_putative_new_pts_3D[t_i].get();

      LandmarkObservations & obs = ml->obs_;

      Vec3 pt_3D_frame_i;
      for(LandmarkObservations::iterator iter_mo = obs.begin(); iter_mo != obs.end();)
      {
        MapObservation & m_o =  iter_mo->second;
        Frame * & frame_i = m_o.frame_ptr;
        // Put 3D point into coordinate system of frame
        getRelativePointPosition(ml->X_,ml->ref_frame_,pt_3D_frame_i,frame_i);

        if (frame_i->getSquaredReprojectionError(pt_3D_frame_i,m_o.feat_id) > frame_i->AC_reprojection_thresh_)
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
        // Delete point from the vector of new points
        #pragma omp critical
        {
          vec_new_pts_outliers.push_back(t_i);
        }
      }
    }
    // sort indexes of outliers
    std::sort(vec_new_pts_outliers.begin(),vec_new_pts_outliers.end());
    // Remove any triangulated landmarks that dont have enough measurements
    for (size_t o_i = 0; o_i < vec_new_pts_outliers.size(); ++o_i)
    {
      std::vector<std::unique_ptr<MapLandmark> >::iterator it_outlier = vec_putative_new_pts_3D.begin();
      // Reduce the number by the number of elements already deleted
      std::advance(it_outlier,vec_new_pts_outliers[o_i]-o_i);

      // Delete 3D point (all other references have been deleted before)
      vec_putative_new_pts_3D.erase(it_outlier);
    }
  }


  // Detect and describe points in current frame
  bool Tracker_Features::detect
  (
    const image::Image<unsigned char> & ima,
    Frame * frame,
    const size_t & min_count,
    const size_t & max_count
  )
  {
    // Detect feature points
    size_t n_feats_detected = featureExtractor_->detect(ima,frame,min_count,max_count);

    if (!(n_feats_detected > 0))
    {
      return false;
    }

    // Describe detected features
    featureExtractor_->describe(ima,frame);

    // Undistort points
    frame->updateFeaturesData();


    if (n_feats_detected > 0)
      return true;
    return false;
  }

} // namespace VO
} // namespace openMVG
