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
    std::shared_ptr<Frame> current_frame,
    const image::Image<unsigned char> * mask
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
    detect(ima,mCurrentFrame.get(),max_tracked_points,mask);

    const size_t & n_feats_detected = mCurrentFrame->getNumberOfFeatures();

    vec_times[1] = omp_get_wtime() - vec_times[0];
    std::cout<<"Tracker: [Detection] #features: "<<n_feats_detected<<" ("<<vec_times[1]<<" s)\n";


    vec_times[2] = omp_get_wtime();
    double start_time = omp_get_wtime();
    // Check if enough features are detected
    if ( ((trackingStatus == TRACKING_STATUS::NOT_INIT || trackingStatus == TRACKING_STATUS::INIT) && n_feats_detected < init_track_min_matches)
        || (trackingStatus == TRACKING_STATUS::OK && n_feats_detected < track_min_matches) )
    {
      // No features detected on the frame
      std::cout<<"Tracker: [Detection] Error - insufficient features detected -> next frame!\n";
      // Reset current frame (keep mPrevFrame from before)
      mCurrentFrame.reset();
      std::cout<<"Tracker: [End tracker] ("<<omp_get_wtime() - start_time<<" s)\n";

      return true;
    }

    // -------------------------
    // -- Tracking
    // -------------------------
    if (trackingStatus == TRACKING_STATUS::NOT_INIT || trackingStatus == TRACKING_STATUS::INIT)
    {
      track_status = tryTrackingInitialization(ima);
    }
    else
    {
      // Flag for marking successful tracking
      bool b_track_OK = false;
      if (trackingStatus == TRACKING_STATUS::OK)
      {
        start_time = omp_get_wtime();
        // ----------------
        // -- Motion Model tracking
        // ----------------
        // Use previous tracking to predict where the pose will be and try matching with previously triangulated map points
        if (motionModel.isValid())
        {
          // Match current frame with lastFrame
          b_track_OK = trackWithMotionModel();
          std::cout<<"Tracker: [Track] Motion Model: "<<b_track_OK<<" ("<<omp_get_wtime() - start_time<<" s)\n";

          display_iterations[0] = 2;
        }

        // If motion model tracking didnt work we try with reference frame
        start_time = omp_get_wtime();
        // ----------------
        // -- Reference frame tracking
        // ----------------
        if (!b_track_OK)
        {
          b_track_OK = trackWithReferenceFrame();
          std::cout<<"Tracker: [Track] Reference frame: "<<b_track_OK<<" ("<<omp_get_wtime() - start_time<<" s)\n";
          display_iterations[1] = 2;
        }
        else
        {
          display_iterations[1] = 0;
        }
        display_iterations[2] = 0;
      }
      else if (trackingStatus == TRACKING_STATUS::LOST)
      {
        start_time = omp_get_wtime();
        // ----------------
        // -- RELOCALIZATION
        // ----------------
        // Try tracking with local map of the last reference frame
        b_track_OK = trackWithReferenceFrameLocalMap(mLastRefFrame.get());

        last_relocalized_frame_id = mCurrentFrame->getFrameId();

        std::cout<<"Tracker: [Track] Relocalization: "<<b_track_OK<<" ("<<omp_get_wtime() - start_time<<" s)\n";

        display_iterations[2] = 2;
      }

      // If successfuly tracked
      if (b_track_OK)
      {
        start_time = omp_get_wtime();
        // ----------------
        // -- Local map tracking
        // ----------------
        trackLocalMap();
        display_iterations[3] = 2;

        size_t n_matches = mCurrentFrame->getNumberMapPoints();
        std::cout<<"Tracker: [Track] Total matches with Map: "<<n_matches<<" ("<<omp_get_wtime() - start_time<<" s)\n";


        mCurrentFrame->computeSceneStatistics();

        start_time = omp_get_wtime();
        // ----------------
        // -- Find new points that we can triangulate
        // ----------------
        std::vector<std::unique_ptr<MapLandmark> > vec_new_pts_3D;
        findNewLandmarks(mCurrentFrame.get(),vec_new_pts_3D);

        std::cout<<"Tracker: [Track Triangulate] #new putative triangulated Pts: "<<vec_new_pts_3D.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";

        ////////////////////////////////////
        // Display new triangulated points
        ///////////////////////////////////
        {
          std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
          for ( size_t k_i = 0; k_i < vec_new_pts_3D.size(); ++k_i)
          {
            std::unique_ptr<MapLandmark> & landmark = vec_new_pts_3D[k_i];

            Vec3 pt_3D_frame;
            PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,mCurrentFrame.get());

            IntrinsicBase * init_cam_intrinsic = mCurrentFrame->getCameraIntrinsics();
            Vec2 pt_3D_frame_projected_init = init_cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());
            d_A.push_back(pt_3D_frame_projected_init);
            d_s.push_back(landmark->getFeatureMeanScale());
            d_C.push_back(Vec2(landmark->association_type_,0));

          }
          display_pt2d_A.push_back(d_A);
          display_pt2d_B.push_back(d_B);
          display_pt2d_C.push_back(d_C);
          display_size_A.push_back(d_s);
          display_text.push_back("Triangulation: Initial new 3D matches");
        }
        ////////////////////////////////////


        start_time = omp_get_wtime();
        if (!cartographer_->optimizeLocal(mCurrentFrame.get(),vec_new_pts_3D,true))
        {
          std::cout<<"Tracker: [Track] Optimize Pose Failed\n";
          return false;
        }

        std::cout<<"Tracker: [Track] Optimize pose OK ("<<omp_get_wtime() - start_time<<" s)\n";

        start_time = omp_get_wtime();
        // ----------------
        // -- Check local landmarks
        // --  - check the reprojection errors and if they have at least 2 inliers
        // ----------------
        cartographer_->verifyLocalLandmarks(mCurrentFrame.get());

        //removeOutliersInLocalFrames(local_map_frames);
        std::cout<<"Tracker: [Track] Local Map verification ("<<omp_get_wtime() - start_time<<" s)\n";

        removeOutliersInFrame(mCurrentFrame.get());
        removeOutliersInNewTriangulatedPoints(mCurrentFrame.get(),vec_new_pts_3D);

        std::cout<<"Tracker: [Track Triangulate] Outlier removal: #new triangulated Pts: "<<vec_new_pts_3D.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";


        ////////////////////////////////////
        // Display new triangulated points
        ///////////////////////////////////
        {
          std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
          for ( size_t k_i = 0; k_i < vec_new_pts_3D.size(); ++k_i)
          {
            std::unique_ptr<MapLandmark> & landmark = vec_new_pts_3D[k_i];

            Vec3 pt_3D_frame;
            PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,mCurrentFrame.get());

            IntrinsicBase * init_cam_intrinsic = mCurrentFrame->getCameraIntrinsics();
            Vec2 pt_3D_frame_projected_init = init_cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());
            d_A.push_back(pt_3D_frame_projected_init);
            d_s.push_back(landmark->getFeatureMeanScale());
            d_C.push_back(Vec2(landmark->association_type_,0));

          }
          display_pt2d_A.push_back(d_A);
          display_pt2d_B.push_back(d_B);
          display_pt2d_C.push_back(d_C);
          display_size_A.push_back(d_s);
          display_text.push_back("Triangulation: Accepted new 3D landmarks");
        }
        ////////////////////////////////////
        display_iterations[5] = 2;



        size_t n_pts_frame = mCurrentFrame->getNumberMapPoints();

        if (n_pts_frame < track_min_matches)
        {
          b_track_OK = false;
        }

        if (b_track_OK)
        {
          // ----------------
          // -- Decide if current frame is keyframe and if should be added to the system
          // ----------------

          bool bKeyframe = needNewKeyframe(mCurrentFrame.get(), vec_new_pts_3D);


          // Add to the system
          if (bKeyframe)
          {
            if (!cartographer_->isMapInitialized())
            {
              // Enough points to add to initialization map
              if(!cartographer_->initializationAddStep(mCurrentFrame, &vec_new_pts_3D))
              {
                std::cout<<"Tracker: [Track] Map Initialization FAILED! Restart everything ("<<omp_get_wtime() - start_time<<" s)\n";
                // Map initialization failed -> reset tracker and map
                cartographer_->clearAllMapData();
                // Reset tracker
                resetTrackingInitialization();
              }
            }
            else
            {
              cartographer_->addStep(mCurrentFrame, &vec_new_pts_3D);
            }

            // Export step to Ply
            std::ostringstream os;
            os << std::setw(8) << std::setfill('0') << "Scene_"<<cartographer_->step_id;
            cartographer_->exportSceneToPly(stlplus::create_filespec("/home/klemen/exp_ply", os.str(), ".ply"), sfm::ESfM_Data((sfm::STRUCTURE | sfm::EXTRINSICS)));

            // Set as reference frame
            mLastRefFrame = mCurrentFrame->share_ptr();
          }

          motionModel.updateMotionModel(mPrevFrame.get(),mCurrentFrame.get());
          track_status = true;
        }
      }
      else
      {
        display_iterations[3] = 0;  // tracking map
        display_iterations[4] = 0;  // epipolar map
        display_iterations[5] = 2;  // triangulations map

      }

      if (!b_track_OK)
      {
        motionModel.setInvalid();
        track_status = false;
        trackingStatus = Abstract_Tracker::TRACKING_STATUS::LOST;

        if (cartographer_->getNumberOfKeyframes() < 5)
        {
          std::cout<<"Tracker: [Track] Got lost at the beginning! Restart everything ("<<omp_get_wtime() - start_time<<" s)\n";
          // Map initialization failed -> reset tracker and map
          cartographer_->clearAllMapData();
          // Reset tracker
          resetTrackingInitialization();
          b_track_OK = false;
          track_status = false;
        }
      }

    }
    // Set current as previous frame
    mPrevPrevFrame.swap(mPrevFrame);
    addPrevFrame(mCurrentFrame);
    mPrevFrame.swap(mCurrentFrame);
    mCurrentFrame.reset();
    std::cout<<"Tracker: [End tracker] ("<<omp_get_wtime() - vec_times[2]<<" s)\n";
    vec_times[3] = omp_get_wtime() - vec_times[2];

    // Return if tracking is ok
    return track_status;
  }

  bool Tracker_Features::needNewKeyframe(Frame * frame,
      std::vector<std::unique_ptr<MapLandmark> > & vec_new_pts_3D)
  {
    // We create new keyframe if:
    // - more than X frames have passed from last mapping
    // - OR
    // - tracking is weak (less than 90% of matches than in last reference frame)
    // - OR
    // - its necessary for creating new points (more than 20% of points are new)

    // Get number of matches with the map
    size_t last_ref_frame_n_matches = mLastRefFrame->getNumberGlobalMapPoints();
    size_t frame_n_matches = frame->getNumberGlobalMapPoints();

    // Check how many frames would become ok if we add this frame as keyframe
    size_t n_new_pts = vec_new_pts_3D.size();

    if (frame_n_matches < last_ref_frame_n_matches*0.9 || n_new_pts > frame_n_matches*0.2 )
      return true;
    else
      return false;
  }

  // --------------------------
  //   Tracking Initialization
  // --------------------------

  // Reset all the variables related to system initialization
  void Tracker_Features::resetTrackingInitialization()
  {
    std::cout<<"Tracker: [Initialization] Reset system initialization process!\n";
    clearTrackingInitializationData();
    std::cout<<"Tracker: [Initialization] Set system status to NOT_INIT\n";
    trackingStatus = Abstract_Tracker::TRACKING_STATUS::NOT_INIT;
    motionModel.setInvalid();

  }

  // Start system initialization process with frame
  void Tracker_Features::startTrackingInitialization(std::shared_ptr<Frame> & frame)
  {
    // Clear all initialization settings
    clearTrackingInitializationData();
    // Set current frame as the new reference frame for initialization
    std::cout<<"Tracker: [Initialization] Set system initialization reference frame\n";
    init_ref_frame = frame->share_ptr();
    // Set system status to INIT
    std::cout<<"Tracker: [Initialization] Set system status to INIT\n";
    trackingStatus = Abstract_Tracker::TRACKING_STATUS::INIT;
  }

  // Try to initialize tracking system
  // Return true if tracking from prev_frame was successful (not necessary initialization)
  bool Tracker_Features::tryTrackingInitialization
  (
    const image::Image<unsigned char> & ima
  )
  {
    if (trackingStatus == Abstract_Tracker::TRACKING_STATUS::NOT_INIT)
    {
      // Set new reference initialization frame
      startTrackingInitialization(mCurrentFrame);
      return false;
    }
    else if (trackingStatus == Abstract_Tracker::TRACKING_STATUS::INIT)
    {
      // Check if we have initialization reference frame
      if (!init_ref_frame)
      {
        // Shouldnt happen
        // If it does we reset the system and and start with this frame as new reference frame
        std::cout<<"Tracker: [Initialization] Missing reference frame! Start initialization with current\n";
        startTrackingInitialization(mCurrentFrame);
        return false;
      }


      double start_time = omp_get_wtime();
      // -------------------
      // -- Match features of two images
      // -------------------
      matching::IndMatches vec_putative_matches_ref_cur_idx;

      // Match all-all features
      featureMatcher_->matching_AllAll_2D_2D(featureExtractor_, init_ref_frame.get(), mCurrentFrame.get(), vec_putative_matches_ref_cur_idx, init_match_desc_ratio,init_match_max_scale_ratio,featureExtractor_->max_dist_desc_);

      std::cout<<"Tracker: [Initialization] [Matching All-All 2D-2D]: " << vec_putative_matches_ref_cur_idx.size() << " (" << omp_get_wtime() - start_time << " s)\n";


      ////////////////////////////////////
      // Display initial pair matches
      ///////////////////////////////////
      {
        std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
        for (auto mp : vec_putative_matches_ref_cur_idx)
        {
          d_A.push_back(init_ref_frame->pts_undist_[mp.i_]);
          d_B.push_back(mCurrentFrame->pts_undist_[mp.j_]);
          d_C.push_back(Vec2(init_ref_frame->pts_scale_[mp.i_],mCurrentFrame->pts_scale_[mp.j_]));
        }
        display_pt2d_A.push_back(d_A);
        display_pt2d_B.push_back(d_B);
        display_pt2d_C.push_back(d_C);
        display_size_A.push_back(d_s);
        display_text.push_back("Initialization: Initial 2D-2D matches");
      }
      ////////////////////////////////////


      // If we dont have enough matches with reference image we set current frame as new reference image
      if (vec_putative_matches_ref_cur_idx.size() < init_track_min_matches)
      {
        std::cout<<"Tracker: [Initialization] Insufficient # matches! Setting this frame as new initialization reference frame\n";
        startTrackingInitialization(mCurrentFrame);
        return false;
      }

      start_time = omp_get_wtime();
      // -------------------
      // -- Estimate relative pose between images
      // -------------------
      bool b_cam_init_calibrated = init_ref_frame->getCamCalibrated();
      bool b_cam_cur_calibrated = mCurrentFrame->getCamCalibrated();

      // Flag for successful estimation of the motion
      bool b_estimated_motion = false;

      // Recovered motion
      Mat4 T_cur_init; // current_T_init
      double AC_reproj_thresh_sq = 16;  // Reprojection error estimated for the model and maximum error
      std::vector<size_t> inliers_M;  // Inliers of an estimated model

      // If calibrated camera we use Essential matrix otherwise Fundamental
      if (b_cam_init_calibrated && b_cam_cur_calibrated)
      {
        b_estimated_motion = PoseEstimator::estimateRobustRelativePosePinholeHE(init_ref_frame.get(), mCurrentFrame.get(), vec_putative_matches_ref_cur_idx, init_min_cos_angle_pt, init_track_min_matches, T_cur_init, inliers_M, AC_reproj_thresh_sq);
      }
      else
      {
        // TODO: Compute HF
        std::cout<<"Tracker: [Initialization] No support for HF initializaton ("<<omp_get_wtime() - start_time<<" s)\n";
      }

      if(!b_estimated_motion)
      {
        std::cout<<"Tracker: [Initialization] Motion estimation failed! Try with next frame\n";
        return true;
      }

      std::cout<<"Tracker: [Initialization] Motion estimation OK: #inliers: "<<inliers_M.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";


      ////////////////////////////////////
      // Display essential matrix inliers
      ///////////////////////////////////
      {
        std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
        for ( size_t k_i = 0; k_i < inliers_M.size(); ++k_i)
        {
          // Advance to the correct inlier
          matching::IndMatches::iterator m_iter = vec_putative_matches_ref_cur_idx.begin();
          std::advance(m_iter,inliers_M[k_i]);
          // Inlier feat ids
          d_A.push_back(init_ref_frame->pts_undist_[m_iter->i_]);
          d_B.push_back(mCurrentFrame->pts_undist_[m_iter->j_]);
          d_C.push_back(Vec2(init_ref_frame->pts_scale_[m_iter->i_],mCurrentFrame->pts_scale_[m_iter->j_]));
        }
        display_pt2d_A.push_back(d_A);
        display_pt2d_B.push_back(d_B);
        display_pt2d_C.push_back(d_C);
        display_size_A.push_back(d_s);
        display_text.push_back("Initialization: Matches with essential matrix");
      }
      ////////////////////////////////////


      // -------------------
      // -- Use estimated pose as initial pose estimation
      // -------------------
      // initial frame
      init_ref_frame->setReferenceFrameAndPose_T(nullptr,Mat4::Identity());
      init_ref_frame->reproj_thresh_sq_ = std::max<double>(default_reproj_thresh_sq,AC_reproj_thresh_sq);
      // current frame
      mCurrentFrame->setReferenceFrameAndPose_T(nullptr,T_cur_init);
      mCurrentFrame->reproj_thresh_sq_ = std::max<double>(default_reproj_thresh_sq,AC_reproj_thresh_sq);

      std::cout<<"Tracker: [Initialization] F_1 thresh: "<<init_ref_frame->reproj_thresh_sq_ <<"\n";
      std::cout<<"Tracker: [Initialization] F_2 thresh: "<<mCurrentFrame->reproj_thresh_sq_ <<"\n";


      start_time = omp_get_wtime();
      // -------------------
      // -- Triangulate inliers
      // -------------------
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
        std::unique_ptr<MapLandmark> & map_landmark = vec_new_triangulated_pts[k_i];
        map_landmark.reset(new MapLandmark());

        // Position of detected features
        //   If we calibrated camera we just use undistorted
        //   If not calibrated we remove the distortion as best as we know it (current estimation)
        const Vec2
          & pt_1 = b_cam_init_calibrated ? init_ref_frame->getFeaturePosition(feat_id_ref) : cam_intrinsic_ref->remove_disto(init_ref_frame->getFeaturePosition(feat_id_ref)),
          & pt_2 = b_cam_cur_calibrated ? mCurrentFrame->getFeaturePosition(feat_id_cur) : cam_intrinsic_cur->remove_disto(mCurrentFrame->getFeaturePosition(feat_id_cur));

        // Triangulate results
        TriangulateDLT(P1, pt_1, P2, pt_2, &(map_landmark->X_));

        // Add observations
        map_landmark->addObservation(init_ref_frame.get(),feat_id_ref);
        map_landmark->addObservation(mCurrentFrame.get(),feat_id_cur);
        map_landmark->setNumberOfObservations(2);

        // Mark as point used in the motion initialization
        map_landmark->association_type_ = 1;
      }
      std::cout<<"Tracker: [Initialization] triangulation ("<<omp_get_wtime() - start_time<<" s)\n";


      ////////////////////////////////////
      // Display triangulated points
      ///////////////////////////////////
      {
        std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
        for ( size_t k_i = 0; k_i < vec_new_triangulated_pts.size(); ++k_i)
        {
          std::unique_ptr<MapLandmark> & landmark = vec_new_triangulated_pts[k_i];

          LandmarkObservations::iterator it_obs = landmark->obs_.begin();

          Vec3 pt_3D_frame;
          PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,(*it_obs).second.frame_ptr);

          IntrinsicBase * init_cam_intrinsic = (*it_obs).second.frame_ptr->getCameraIntrinsics();
          Vec2 pt_3D_frame_projected_init = init_cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());
          d_A.push_back(pt_3D_frame_projected_init);

          it_obs++;

          PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,(*it_obs).second.frame_ptr);
          init_cam_intrinsic = (*it_obs).second.frame_ptr->getCameraIntrinsics();
          pt_3D_frame_projected_init = init_cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());
          d_B.push_back(pt_3D_frame_projected_init);
          d_C.push_back(Vec2(0,0));
        }
        display_pt2d_A.push_back(d_A);
        display_pt2d_B.push_back(d_B);
        display_pt2d_C.push_back(d_C);
        display_size_A.push_back(d_s);
        display_text.push_back("Initialization: Triangulated points");
      }
      ////////////////////////////////////


      start_time = omp_get_wtime();
      // -------------------
      // -- Optimize initial poses
      // -------------------
      // First frame is fixed
      init_ref_frame->setActive();

      std::cout<<"Tracker: [Initialization] # of initial tracked points before: "<<vec_new_triangulated_pts.size()<<"\n";

      // Optimize first pair
      if (!cartographer_->optimizeLocal(mCurrentFrame.get(), vec_new_triangulated_pts,true))
      {
        std::cout<<"Tracker: [Initialization] Optimize local FAILED ("<<omp_get_wtime() - start_time<<" s)\n";
        return false;
      }

      std::cout<<"Tracker: [Initialization] Optimize local OK ("<<omp_get_wtime() - start_time<<" s)\n";


      start_time = omp_get_wtime();
      // -------------------
      // -- Determine inliers
      // -------------------

      // Compute scene borders
      //init_ref_frame->computeSceneStatistics();
      //mCurrentFrame->computeSceneStatistics();
      // Check for outliers in the frames
      removeOutliersInNewTriangulatedPoints(mCurrentFrame.get(),vec_new_triangulated_pts);


      std::cout<<"Tracker: [Initialization] # of initial tracked points: "<<vec_new_triangulated_pts.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";

      ////////////////////////////////////
      // Display triangulated points - INLIERS
      ///////////////////////////////////
      {
        std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
        for ( size_t k_i = 0; k_i < vec_new_triangulated_pts.size(); ++k_i)
        {
          std::unique_ptr<MapLandmark> & landmark = vec_new_triangulated_pts[k_i];

          LandmarkObservations::iterator it_obs = landmark->obs_.begin();

          Vec3 pt_3D_frame;
          PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,(*it_obs).second.frame_ptr);

          IntrinsicBase * init_cam_intrinsic = (*it_obs).second.frame_ptr->getCameraIntrinsics();
          Vec2 pt_3D_frame_projected_init = init_cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());
          d_A.push_back(pt_3D_frame_projected_init);

          it_obs++;

          PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,(*it_obs).second.frame_ptr);
          init_cam_intrinsic = (*it_obs).second.frame_ptr->getCameraIntrinsics();
          pt_3D_frame_projected_init = init_cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());
          d_B.push_back(pt_3D_frame_projected_init);
          d_C.push_back(Vec2(0,0));
        }
        display_pt2d_A.push_back(d_A);
        display_pt2d_B.push_back(d_B);
        display_pt2d_C.push_back(d_C);
        display_size_A.push_back(d_s);
        display_text.push_back("Initialization: Accepted triangulated points");
      }
      ////////////////////////////////////


      // -------------------
      // -- Check that there are enough tracked points
      // -------------------
      if (vec_new_triangulated_pts.size() < init_track_min_matches)
      {
        // Not enough matches with initialization reference frame
        // Set this frame as initialization reference frame and try with next frame
        std::cout<<"Tracker: [Initialization] Insufficient # of tracked points: "<<vec_new_triangulated_pts.size()<<"! Setting this frame as new initialization reference frame\n";
        startTrackingInitialization(mCurrentFrame);
        return false;
      }

      start_time = omp_get_wtime();
      // -------------------
      // -- Add first two frames to map
      // -------------------
      cartographer_->initializationAddStep(init_ref_frame, nullptr);
      cartographer_->initializationAddStep(mCurrentFrame, &vec_new_triangulated_pts);

      std::cout<<"Tracker: [Initialization] Frames 1,2 added to the map ("<<omp_get_wtime() - start_time<<" s)\n";

      // Update motion model
      motionModel.updateMotionModel(init_ref_frame.get(),mCurrentFrame.get());

      // Set second camera as reference keyframe
      mLastRefFrame = mCurrentFrame->share_ptr();

      // Clear initialization data
      clearTrackingInitializationData();

      // Set tracking to OK
      std::cout<<"Tracker: [Initialization] Set system status to OK\n";trackingStatus = Abstract_Tracker::TRACKING_STATUS::OK;
      return true;
    }
  }


  // --------------------------
  //   Tracking
  // --------------------------
  bool Tracker_Features::trackWithMotionModel()
  {
    Hash_Map<MapLandmark *,IndexT> map_putative_matches_3D_pts_frame_idx;

    std::cout<<"Tracker: [Track] Use Motion Model\n";

    // -------------------
    // -- Predict location of current frame (using MotionModel)
    // -------------------
    Mat4 T_predict = motionModel.predict(mPrevFrame.get(),mCurrentFrame.get());
    // TODO: Check how to handle relative poses
    // As motion model has forward transformations (we have to inverte before setting it)
    mCurrentFrame->setPose_T(T_predict.inverse(), nullptr);

    double start_time = omp_get_wtime();
    // -------------------
    // -- Match by projecting triangulated points from last track_mm_n_frames frames to current frame
    // -------------------
    featureMatcher_->matching_Projection_3D_2D(featureExtractor_, mPrevFrame.get(),mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx, track_mm_win_size, track_match_mm_desc_ratio, track_match_mm_max_scale_ratio, featureExtractor_->max_dist_desc_);

    std::cout<<"Tracker: [Track] Match by projection with frame: "<<mPrevFrame->getFrameId()<<": Total # matches: "<<map_putative_matches_3D_pts_frame_idx.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";


    // If not enough matches we try again with wider search window
    if (map_putative_matches_3D_pts_frame_idx.size() < track_min_matches)
    {
      start_time = omp_get_wtime();
      // Double search window for searching by projection
      featureMatcher_->matching_Projection_3D_2D(featureExtractor_, mPrevFrame.get(),mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx, track_mm_win_size*2, track_match_mm_desc_ratio, track_match_mm_max_scale_ratio, featureExtractor_->max_dist_desc_);

      std::cout<<"Tracker: [Track] Match by projection with frame: "<<mPrevFrame->getFrameId()<<": Total # matches: "<<map_putative_matches_3D_pts_frame_idx.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";


      if (map_putative_matches_3D_pts_frame_idx.size() < track_min_matches)
      {
        std::cout<<"Tracker: [Track] Match by Motion Model Failed\n";
        return false;
      }
    }

    ////////////////////////////////////
    // Display initial matches proj 3d-2d using motion model
    ///////////////////////////////////
    {
      IntrinsicBase * cam_intrinsic = mCurrentFrame->getCameraIntrinsics();
      std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
      for ( size_t k_i = 0; k_i < mPrevFrame->map_points_.size(); ++k_i)
      {
        MapLandmark * landmark = mPrevFrame->map_points_[k_i];

        if (!landmark)
          continue;

        Vec3 pt_3D_frame;
        PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,mCurrentFrame.get());

        Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());

        d_B.push_back(mPrevFrame->pts_undist_[k_i]);

        if (map_putative_matches_3D_pts_frame_idx.find(landmark)!=map_putative_matches_3D_pts_frame_idx.end())
        {
          d_C.push_back(mCurrentFrame->pts_undist_[map_putative_matches_3D_pts_frame_idx.find(landmark)->second]);
        }
        else
        {
          d_C.push_back(Vec2(-1,-1));
        }

        d_A.push_back(pt_3D_frame_projected);

        // Compute viewing angle between the point and normal which was last seen
        Vec3 normal_X = (landmark->getWorldPosition() - mCurrentFrame->O_w_).normalized();
        const Vec3 & normal_last = landmark->getLastNormal();

        // Compute viewing angle
        float angle_factor = featureMatcher_->radiusByViewingAngle(normal_X.dot(normal_last));   // [1,5]

        // enlarge the area if the viewing angle is bigger
        float win_size = angle_factor;//track_mm_win_size * sqrt(angle_factor);

        d_s.push_back(win_size);
      }
      display_pt2d_A.push_back(d_A);
      display_pt2d_B.push_back(d_B);
      display_pt2d_C.push_back(d_C);
      display_size_A.push_back(d_s);
      display_text.push_back("MotionModel: Initial 3D-2D matches");
    }
    ////////////////////////////////////


    // -------------------
    // -- Do tiny BA to refine the pose based on the matches
    // -------------------
    mCurrentFrame->reproj_thresh_sq_ = default_reproj_thresh_sq;

    start_time = omp_get_wtime();
    // Optimize pose with static points (both local and global)
    if (!cartographer_->optimizePose(mCurrentFrame.get(),map_putative_matches_3D_pts_frame_idx,true))
    {
      std::cout<<"Tracker: [Track] Optimize Pose Failed\n";
      return false;
    }

    std::cout<<"Tracker: [Track] Optimize pose OK ("<<omp_get_wtime() - start_time<<" s)\n";

    removeOutliersInFrame(mCurrentFrame.get());
    // Check Chi2 error for each point and mark which ones are ok
    markInliersInFrame(mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx,true,2);
    size_t n_matches = mCurrentFrame->getNumberMapPoints();

    std::cout<<"Tracker: [Track] # matched with motion model: "<<n_matches<<" ("<<omp_get_wtime() - start_time<<" s)\n";


    ////////////////////////////////////
    // Display matches proj 3d-2d using motion model after BA
    ///////////////////////////////////
    {
      IntrinsicBase * cam_intrinsic = mCurrentFrame->getCameraIntrinsics();
      std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
      for ( size_t k_i = 0; k_i < mPrevFrame->map_points_.size(); ++k_i)
      {
        MapLandmark * landmark = mPrevFrame->map_points_[k_i];

        if (!landmark)
          continue;

        Vec3 pt_3D_frame;
        PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,mCurrentFrame.get());

        Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());

        d_B.push_back(mPrevFrame->pts_undist_[k_i]);

        if (map_putative_matches_3D_pts_frame_idx.find(landmark)!=map_putative_matches_3D_pts_frame_idx.end())
        {
          d_C.push_back(mCurrentFrame->pts_undist_[map_putative_matches_3D_pts_frame_idx.find(landmark)->second]);
        }
        else
        {
          d_C.push_back(Vec2(-1,-1));
        }

        d_A.push_back(pt_3D_frame_projected);

        // Compute viewing angle between the point and normal which was last seen
        Vec3 normal_X = (landmark->getWorldPosition() - mCurrentFrame->O_w_).normalized();
        const Vec3 & normal_last = landmark->getLastNormal();

        // Compute viewing angle
        float angle_factor = featureMatcher_->radiusByViewingAngle(normal_X.dot(normal_last));   // [1,5]

        // enlarge the area if the viewing angle is bigger
        float win_size = track_mm_win_size * sqrt(angle_factor);

        d_s.push_back(win_size);
      }
      display_pt2d_A.push_back(d_A);
      display_pt2d_B.push_back(d_B);
      display_pt2d_C.push_back(d_C);
      display_size_A.push_back(d_s);
      display_text.push_back("MotionModel: Accepted 3D-2D matches after BA");
    }
    ////////////////////////////////////


    // -------------------
    // -- Check if motion model matching is successful
    // -------------------
    if (n_matches < track_min_matches)
    {
      return false;
    }
    return true;
  }

  bool Tracker_Features::trackWithReferenceFrame()
  {
    std::cout<<"Tracker: [Track] Use Reference frame\n";
    // -------------------
    // -- Set location of current frame (as the one of last frame)
    // -------------------
    // TODO: check how handle relative poses
    mCurrentFrame->setPose_T(mPrevFrame->getTransformationMatrix(), nullptr);

    double start_time = omp_get_wtime();
    // -------------------
    // -- Try to match all-all with the last reference frame
    // -------------------
    Hash_Map<MapLandmark* ,IndexT> map_putative_matches_3D_pts_frame_idx;
    featureMatcher_->matching_AllAll_3D_2D(featureExtractor_, mLastRefFrame->map_points_, mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx, track_match_rf_desc_ratio, track_match_rf_max_scale_ratio, featureExtractor_->max_dist_desc_);

    std::cout<<"Tracker: [Track] All-All 3D-2D! RF ID: "<<mLastRefFrame->getFrameId()<<" with # pts: "<<mLastRefFrame->getNumberMapPoints()<<"\n";
    std::cout<<"Tracker: [Track] All-All 3D-2D! # matches: "<<map_putative_matches_3D_pts_frame_idx.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";



    ////////////////////////////////////
    // Display initial matches proj 3d-2d using motion model
    ///////////////////////////////////
    {
      IntrinsicBase * cam_intrinsic = mCurrentFrame->getCameraIntrinsics();
      std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
      for ( size_t k_i = 0; k_i < mLastRefFrame->map_points_.size(); ++k_i)
      {
        MapLandmark * landmark = mLastRefFrame->map_points_[k_i];

        if (!landmark)
          continue;

        Vec3 pt_3D_frame;
        PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,mCurrentFrame.get());

        Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());

        d_B.push_back(mLastRefFrame->pts_undist_[k_i]);

        if (map_putative_matches_3D_pts_frame_idx.find(landmark)!=map_putative_matches_3D_pts_frame_idx.end())
        {
          d_C.push_back(mCurrentFrame->pts_undist_[map_putative_matches_3D_pts_frame_idx.find(landmark)->second]);
        }
        else
        {
          d_C.push_back(Vec2(-1,-1));
        }

        d_A.push_back(pt_3D_frame_projected);


        // enlarge the area if the viewing angle is bigger
        float win_size = 1;

        d_s.push_back(win_size);
      }
      display_pt2d_A.push_back(d_A);
      display_pt2d_B.push_back(d_B);
      display_pt2d_C.push_back(d_C);
      display_size_A.push_back(d_s);
      display_text.push_back("ReferenceFrame: Inital 3D-2D matches");
    }
    ////////////////////////////////////


    // Check if we got enough matches
    if (map_putative_matches_3D_pts_frame_idx.size() < track_min_matches)
    {
      std::cout<<"Tracker: [Track] Match with Reference Frame Failed\n";
      return false;
    }

    // TODO: Check how to set reprojection threshold if we just use BA
    mCurrentFrame->reproj_thresh_sq_ = default_reproj_thresh_sq;

    // Optimize pose with static points (both local and global)
    start_time = omp_get_wtime();
    if (!cartographer_->optimizePose(mCurrentFrame.get(),map_putative_matches_3D_pts_frame_idx,true))
    {
      std::cout<<"Tracker: [Track] Optimize Pose Failed\n";
      return false;
    }

    std::cout<<"Tracker: [Track] Optimize pose OK ("<<omp_get_wtime() - start_time<<" s)\n";

    removeOutliersInFrame(mCurrentFrame.get());
    markInliersInFrame(mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx,true,3);

    size_t n_matches = mCurrentFrame->getNumberMapPoints();

    std::cout<<"Tracker: [Track] # matched with reference frame: "<<n_matches<<" ("<<omp_get_wtime() - start_time<<" s)\n";


    ////////////////////////////////////
    // Display link between projeccted point and feature on the frame
    ///////////////////////////////////
    {
      IntrinsicBase * cam_intrinsic = mCurrentFrame->getCameraIntrinsics();

      std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
      for ( size_t k_i = 0; k_i < mCurrentFrame->map_points_.size(); ++k_i)
      {
        MapLandmark * landmark = mCurrentFrame->map_points_[k_i];

        if (!landmark)
          continue;

        Vec3 pt_3D_frame;
        PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,mCurrentFrame.get());

        Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());
        d_B.push_back(Vec2(-1,-1));
        d_C.push_back(mCurrentFrame->pts_undist_[k_i]);
        d_A.push_back(pt_3D_frame_projected);

        Hash_Map<MapLandmark *,IndexT>::iterator it = map_putative_matches_3D_pts_frame_idx.find(landmark);
        size_t it_dist = std::distance(map_putative_matches_3D_pts_frame_idx.begin(), it);

        d_s.push_back(display_size_A[display_size_A.size()-1][it_dist]);
      }
      display_pt2d_A.push_back(d_A);
      display_pt2d_B.push_back(d_B);
      display_pt2d_C.push_back(d_C);
      display_size_A.push_back(d_s);
      display_text.push_back("ReferenceFrame: Accepted 3D-2D matches after BA");
    }
    ////////////////////////////////////

    // -------------------
    // -- Check if motion model matching is successful
    // -------------------
    if (n_matches < track_min_matches)
    {
      return false;
    }
    return true;



  }


  bool Tracker_Features::trackWithReferenceFrameLocalMap(Frame * frame_ref)
  {
    std::cout<<"Tracker: [Track] Use Local Map of Reference Frame\n";

    double start_time = omp_get_wtime();
    // -------------------
    // -- Try to match 3D points of a local map of last reference frame with current frame
    // -------------------
    // -- Identify local map (frames & points)
    std::vector<Frame *> local_map_frames;
    std::vector<MapLandmark *> local_map_points;

    // Identify [N = track_local_map_size] best frames that see the most points that are also seen from current camera
    frame_ref->getFrameVisibilityConnections(local_map_frames,10);

    // Get all unmatched map points from frames in local map
    cartographer_->getLocalMapPoints(mCurrentFrame.get(),local_map_frames,local_map_points);

    // Set local point back to not used
    cartographer_->resetFlagLocalMapPoints(local_map_points);

    std::cout<<"Tracker: [Track] Get Frame Visibility. # local frames: " << local_map_frames.size()<<" # local pts: "<<local_map_points.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";


    start_time = omp_get_wtime();
    // -------------------
    // -- Try to match local map points with the features in the current frame
    // -------------------
    Hash_Map<MapLandmark *,IndexT> map_putative_matches_3D_pts_frame_idx_putative;


    // Matching All 3D points in the local map with all 2D features in the image
    featureMatcher_->matching_AllAll_3D_2D(featureExtractor_, local_map_points, mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx_putative, track_match_rm_desc_ratio,track_match_rm_max_scale_ratio, featureExtractor_->max_dist_desc_);

    std::cout<<"Tracker: [Track] Match All 3D - 2D  #matches: "<< map_putative_matches_3D_pts_frame_idx_putative.size()<<"\n";

    // Estimate the pose based on the 3D - 2D matching (use scale from last reference frame
    Mat3 R;Vec3 t; double s=frame_ref->getTransformationMatrix().block(0,0,3,1).norm();
    std::vector<size_t> vec_inliers;
    double AC_threshold;

    bool b_resection = PoseEstimator::estimateRobustPose(mCurrentFrame.get(),map_putative_matches_3D_pts_frame_idx_putative,track_min_matches,R,t,vec_inliers,AC_threshold);

    // Check if we got enough matches
    if (!b_resection || vec_inliers.size() < track_min_matches)
    {
      std::cout<<"Tracker: [Track] Match with LocalMap Failed\n";
      return false;
    }

    // Save inliers
    Hash_Map<MapLandmark *,IndexT> map_putative_matches_3D_pts_frame_idx;
    for (auto inlier_i: vec_inliers)
    {
      Hash_Map<MapLandmark *,IndexT>::iterator it = map_putative_matches_3D_pts_frame_idx_putative.begin();
      std::advance(it,inlier_i);
      map_putative_matches_3D_pts_frame_idx[it->first] = it->second;
    }
    std::cout<<"Tracker: [Track] Matches with Local Map of Reference frame #matches: "<<map_putative_matches_3D_pts_frame_idx.size()<<"\n";

    // Set pose
    mCurrentFrame->setPose_Rts(R,t,s,nullptr);

    ////////////////////////////////////
    // Display initial matches proj 3d-2d using motion model
    ///////////////////////////////////
    {
      IntrinsicBase * cam_intrinsic = mCurrentFrame->getCameraIntrinsics();
      std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
      for ( size_t k_i = 0; k_i < local_map_points.size(); ++k_i)
      {
        MapLandmark * landmark = local_map_points[k_i];

        if (!landmark)
          continue;

        Vec3 pt_3D_frame;
        PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,mCurrentFrame.get());

        Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());

        if (!mCurrentFrame->isPointInFrame(pt_3D_frame_projected))
        {
          continue;
        }

        d_B.push_back(Vec2(-1,-1));

        if (map_putative_matches_3D_pts_frame_idx.find(landmark)!=map_putative_matches_3D_pts_frame_idx.end())
        {
          d_C.push_back(mCurrentFrame->pts_undist_[map_putative_matches_3D_pts_frame_idx.find(landmark)->second]);
        }
        else
        {
          d_C.push_back(Vec2(-1,-1));
        }

        d_A.push_back(pt_3D_frame_projected);


        // enlarge the area if the viewing angle is bigger
        float win_size = 10;

        d_s.push_back(win_size);
      }
      display_pt2d_A.push_back(d_A);
      display_pt2d_B.push_back(d_B);
      display_pt2d_C.push_back(d_C);
      display_size_A.push_back(d_s);
      display_text.push_back("ReferenceMap: Initial 3D-2D matches");
    }
    ////////////////////////////////////


    // TODO: Check how to set reprojection threshold if we just use BA
    mCurrentFrame->reproj_thresh_sq_ = default_reproj_thresh_sq;

    // Optimize pose with static points (both local and global)
    start_time = omp_get_wtime();
    if (!cartographer_->optimizePose(mCurrentFrame.get(),map_putative_matches_3D_pts_frame_idx,true))
    {
      std::cout<<"Tracker: [Track] Optimize Pose Failed\n";
      return false;
    }

    std::cout<<"Tracker: [Track] Optimize pose OK ("<<omp_get_wtime() - start_time<<" s)\n";

    removeOutliersInFrame(mCurrentFrame.get());
    markInliersInFrame(mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx,true,4);
    size_t n_matches = mCurrentFrame->getNumberMapPoints();

    std::cout<<"Tracker: [Track] # matched with local map of reference map: "<<n_matches<<" ("<<omp_get_wtime() - start_time<<" s)\n";

    ////////////////////////////////////
    // Display link between projected point and feature on the frame
    ///////////////////////////////////
    {
      IntrinsicBase * cam_intrinsic = mCurrentFrame->getCameraIntrinsics();

      std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
      for ( size_t k_i = 0; k_i < mCurrentFrame->map_points_.size(); ++k_i)
      {
        MapLandmark * landmark = mCurrentFrame->map_points_[k_i];

        if (!landmark)
          continue;

        Vec3 pt_3D_frame;
        PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,mCurrentFrame.get());

        Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());
        d_B.push_back(Vec2(-1,-1));
        d_C.push_back(mCurrentFrame->pts_undist_[k_i]);
        d_A.push_back(pt_3D_frame_projected);

        Hash_Map<MapLandmark *,IndexT>::iterator it = map_putative_matches_3D_pts_frame_idx.find(landmark);
        size_t it_dist = std::distance(map_putative_matches_3D_pts_frame_idx.begin(), it);

        d_s.push_back(display_size_A[display_size_A.size()-1][it_dist]);
      }
      display_pt2d_A.push_back(d_A);
      display_pt2d_B.push_back(d_B);
      display_pt2d_C.push_back(d_C);
      display_size_A.push_back(d_s);
      display_text.push_back("ReferenceMap: Accepted 3D-2D matches after BA");
    }
    ////////////////////////////////////

    // -------------------
    // -- Check if motion model matching is successful
    // -------------------
    if (n_matches < track_min_matches)
    {
      return false;
    }
    return true;
  }


  bool Tracker_Features::trackLocalMap()
  {
    std::cout<<"Tracker: [Track] Local map\n";

    double start_time = omp_get_wtime();
    // -------------------
    // -- Identify local map (frames & points)
    // -------------------
    std::vector<Frame *> local_map_frames;
    std::vector<MapLandmark *> local_map_points;

    // Identify [N = track_local_map_size] best frames that see the most points that are also seen from current camera
    mCurrentFrame->getFrameVisibilityConnections(local_map_frames,track_local_map_size);

    // Get all unmatched map points from frames in local map
    cartographer_->getLocalMapPoints(mCurrentFrame.get(),local_map_frames,local_map_points);

    std::cout<<"Tracker: [Track Local Map] Get Frame Visibility. # local frames: " << local_map_frames.size()<<" # local pts: "<<local_map_points.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";


    start_time = omp_get_wtime();
    // -------------------
    // -- Try to match local map points with the features in the current frame
    // -------------------
    Hash_Map<MapLandmark *,IndexT> map_putative_matches_3D_pts_frame_idx;
    // If we did relocalization recently we double the size of the window
    float t_lm_win_size = track_lm_win_size;
    if (mCurrentFrame->getFrameId() < last_relocalized_frame_id+2)
    {
      t_lm_win_size *= 2;
    }

    featureMatcher_->matching_Projection_3D_2D(featureExtractor_, mCurrentFrame.get(), local_map_points, map_putative_matches_3D_pts_frame_idx, t_lm_win_size, track_match_lm_desc_ratio, track_match_lm_max_scale_ratio, featureExtractor_->max_dist_desc_);

    std::cout<<"Tracker: [Track Local Map] Match by projection 3D-2D with Local Map! # matches: "<<map_putative_matches_3D_pts_frame_idx.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";

    ////////////////////////////////////
    // Display initial matches proj 3d-2d using local map
    ///////////////////////////////////
    {
      std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
      for ( size_t k_i = 0; k_i < local_map_points.size(); ++k_i)
      {
        MapLandmark * landmark = local_map_points[k_i];

        if (!landmark)
          continue;

        Vec3 pt_3D_frame;
        PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,mCurrentFrame.get());

        if (pt_3D_frame(2) < 0)
          continue;

        IntrinsicBase * cam_intrinsic = mCurrentFrame->getCameraIntrinsics();
        Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());
        // Check if projection is actually in the image borders
        if (!mCurrentFrame->isPointInFrame(pt_3D_frame_projected))
          continue;

        d_B.push_back(Vec2(-1,-1));
        if (map_putative_matches_3D_pts_frame_idx.find(landmark)!=map_putative_matches_3D_pts_frame_idx.end())
        {
          d_C.push_back(mCurrentFrame->pts_undist_[map_putative_matches_3D_pts_frame_idx.find(landmark)->second]);
        }
        else
        {
          d_C.push_back(Vec2(-1,-1));
        }

        d_A.push_back(pt_3D_frame_projected);

        d_s.push_back(t_lm_win_size);
      }
      display_pt2d_A.push_back(d_A);
      display_pt2d_B.push_back(d_B);
      display_pt2d_C.push_back(d_C);
      display_size_A.push_back(d_s);
      display_text.push_back("LocalMap: Initial 3D-2D matches");
    }
    ////////////////////////////////////

    start_time = omp_get_wtime();
    // ----------------
    // -- Perform BA on all the matches existing points
    // ----------------
    double lm_inlier_ratio = computeInlierRatioMatches(mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx);

    //if (lm_inlier_ratio < track_lm_pose_BA_inlier_ratio)
    {
      std::cout<<"Tracker: [Track] Inlier ratio of matched map landmarks ( "<<lm_inlier_ratio<<" ) is below threshold! Perform pose BA.\n";

      start_time = omp_get_wtime();
      // Optimize pose with static points (both local and global)
      if (!cartographer_->optimizePose(mCurrentFrame.get(),map_putative_matches_3D_pts_frame_idx,true))
      {
        std::cout<<"Tracker: [Track] Optimize Pose Failed\n";
        return false;
      }

      std::cout<<"Tracker: [Track] Optimize pose OK ("<<omp_get_wtime() - start_time<<" s)\n";
    }

    removeOutliersInFrame(mCurrentFrame.get());
    // Add matches as inliers -> we will check if they are outliers after local optimization
    markInliersInFrame(mCurrentFrame.get(), map_putative_matches_3D_pts_frame_idx,false,5);


    ////////////////////////////////////
    // Display link between projeccted point and feature on the frame
    ///////////////////////////////////
    {
      IntrinsicBase * cam_intrinsic = mCurrentFrame->getCameraIntrinsics();

      std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;
      for ( size_t k_i = 0; k_i < mCurrentFrame->map_points_.size(); ++k_i)
      {
        MapLandmark * landmark = mCurrentFrame->map_points_[k_i];

        if (!landmark)
          continue;

        Vec3 pt_3D_frame;
        PoseEstimator::getRelativePointPosition(landmark->X_,landmark->ref_frame_,pt_3D_frame,mCurrentFrame.get());

        Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());
        d_B.push_back(Vec2(-1,-1));
        d_C.push_back(mCurrentFrame->pts_undist_[k_i]);
        d_A.push_back(pt_3D_frame_projected);

        Hash_Map<MapLandmark *,IndexT>::iterator it = map_putative_matches_3D_pts_frame_idx.find(landmark);
        size_t it_dist = std::distance(map_putative_matches_3D_pts_frame_idx.begin(), it);

        d_s.push_back(display_size_A[display_size_A.size()-1][it_dist]);
      }
      display_pt2d_A.push_back(d_A);
      display_pt2d_B.push_back(d_B);
      display_pt2d_C.push_back(d_C);
      display_size_A.push_back(d_s);
      display_text.push_back("LocalMap: Accepted 3D-2D matches after BA");
    }
    ////////////////////////////////////
    display_iterations.push_back(2);
    return true;
  }

  void Tracker_Features::findNewLandmarks(Frame * frame, std::vector<std::unique_ptr<MapLandmark> > & vec_new_pts_3D)
  {
    std::vector<Frame *> local_map_frames;

    // Compute scene statistics of the focus frame
    frame->computeSceneStatistics();
    // Center of focus camera
    const Vec3 frame_camera_center = frame->getCameraCenter();

    double start_time = omp_get_wtime();
    // -------------------
    // -- Identify local map (frames)
    // -------------------
    // Identify N best frames that see points seen from current camera
    frame->getFrameVisibilityConnections(local_map_frames,triangulate_local_map_size);

    std::cout<<"Tracker: [Track Triangulate] Get Frame Visibility. # local frames: " << local_map_frames.size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";

    // -------------------
    // -- Loop through each neighbor frame and check if we can match any of the (unmatched) points in current frame to
    // -- (unmatched) features in the neighbor frames
    // -------------------
    display_iterations[4] = local_map_frames.size();
    std::vector<Hash_Map<IndexT,IndexT> > vec_matches_cur_local(local_map_frames.size());

    start_time = omp_get_wtime();
    std::vector<float> vec_scene_max(local_map_frames.size(),0.0);
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t f_i = 0; f_i < local_map_frames.size(); ++f_i)
    {
      Frame * local_frame_i = local_map_frames[f_i];

      // Compute ratio baseline/scene depth between the pair of frames
      const Vec3 frame_i_camera_center = local_frame_i->getCameraCenter();
      local_frame_i->computeSceneStatistics();
      // The limit for new points will be 3 times the 95% of the points distance

      const float baseline = (frame_camera_center - frame_i_camera_center).norm();
      const float ratio = baseline / local_frame_i->f_scene_median;
      std::cout<<"Tracker: [Track Triangulate] Ratio baseline/scene depth: "<<ratio<< " (baseline: "<<baseline<<" depth (median): "<<local_frame_i->f_scene_median<<")\n";

      if (ratio > 0.01)
      {
        // Compute fundamental matrix between frame_current and F_local
        Mat3 F_l_c; // p_local' l_F_c p_cur
        PoseEstimator::computeFundamentalMatrix(frame,local_frame_i,F_l_c);

        // Perform epipolar matching
        featureMatcher_->matching_EpipolarLine_2D_2D(featureExtractor_, frame, local_frame_i, F_l_c, vec_matches_cur_local[f_i],track_epipolar_dist,triangulate_match_epipolar_desc_ratio,triangulate_match_lm_max_scale_ratio, featureExtractor_->max_dist_desc_);

        std::cout<<"Tracker: [Track Triangulate] Match by epipolar 2D-2D with frame: "<<local_frame_i->getFrameId()<<"! # matches: "<<vec_matches_cur_local[f_i].size()<<" ("<<omp_get_wtime() - start_time<<" s)\n";

        ////////////////////////////////////
        // Display epipolar mappings
        ///////////////////////////////////
        {
          std::vector<Vec2> d_A; std::vector<Vec2> d_B;std::vector<Vec2> d_C; std::vector<float> d_s;

          Vec3 pt_epi_3D;
         PoseEstimator::getRelativePointPosition(frame_camera_center,nullptr,pt_epi_3D,local_frame_i);

         // Project the point to image plane of frame_2 2
         const Vec2 pt_epi_3D_proj = local_frame_i->getCameraIntrinsics()->cam2ima(local_frame_i->getCameraIntrinsics()->have_disto()?local_frame_i->getCameraIntrinsics()->add_disto(pt_epi_3D.hnormalized()):pt_epi_3D.hnormalized());

          d_A.push_back(pt_epi_3D_proj);
          d_B.push_back(pt_epi_3D_proj);
          d_s.push_back(1);

          for ( auto match : vec_matches_cur_local[f_i])
          {

            d_A.push_back(frame->pts_undist_[match.first]);
            d_B.push_back(local_frame_i->pts_undist_[match.second]);

            d_s.push_back(local_frame_i->pts_scale_[match.second]);
          }
          display_pt2d_A.push_back(d_A);
          display_pt2d_B.push_back(d_B);
          display_pt2d_C.push_back(d_C);
          display_size_A.push_back(d_s);
          display_text.push_back("EpipolarMatch: Initial 2D-2D matches");
        }
        ////////////////////////////////////

      }
      else
      {
        std::cout<<"Tracker: [Track Triangulate] Skip pair. Small ratio baseline/scene depth: "<<ratio<< " (baseline: "<<baseline<<" depth (median): "<<local_frame_i->f_scene_median<<")\n";
      }
    }

    std::cout<<"Tracker: [Track Triangulate] Epipolar 2D 2D: "<<omp_get_wtime() - start_time<<" s\n";

    const size_t frame_n_feats = frame->getNumberOfFeatures();
    // Reserve max possible new points
    vec_new_pts_3D.reserve(frame_n_feats);

    start_time = omp_get_wtime();
    // -------------------
    // -- Concatenate observations of the same point over all neighbor frames and triangulate
    // -------------------
    std::vector<bool> b_feat_searched = std::vector<bool>(frame_n_feats,false);

    // Projection matrix
    Mat34 P,P_i;

    // Information for focus frame
    const Mat3 frame_Rwc = frame->getRotationMatrixInverse();
    const IntrinsicBase * cam_intrinsic = frame->getCameraIntrinsics();
    const Mat3 K_inv = frame->getK_inv();
    //TODO: How to use if we have relative cameras
    P = frame->getProjectionMatrix();


    // Go through each local frame and add concatenate observations per feature in current frame
    for (size_t local_frame_i = 0; local_frame_i < local_map_frames.size(); ++local_frame_i)
    {
      Frame * local_frame = local_map_frames[local_frame_i];
      const IndexT local_frame_id = local_frame->getFrameId();

      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for (size_t m_i = 0; m_i < vec_matches_cur_local[local_frame_i].size(); ++m_i)
      {
        Hash_Map<IndexT,IndexT>::iterator match_i = vec_matches_cur_local[local_frame_i].begin();
        std::advance(match_i,m_i);
        IndexT feat_id = match_i->first;

        // If we already searched for the feature we skip as we have already searched for it
        if (b_feat_searched[feat_id])
          continue;

        // Create the concatenation of observations from other frames
        LandmarkObservations feat_obs_nviews;
        feat_obs_nviews[local_frame_id] = MapObservation(match_i->second, local_frame);

        // We only have to look for next  frames as if it would exist previously we would already process it
        for (size_t local_frame_j = local_frame_i+1; local_frame_j < local_map_frames.size(); ++local_frame_j)
        {
          Hash_Map<IndexT,IndexT>::iterator it_matches_cur_frame_j = vec_matches_cur_local[local_frame_j].find(feat_id);
          if (it_matches_cur_frame_j != vec_matches_cur_local[local_frame_j].end())
          {
            feat_obs_nviews[local_map_frames[local_frame_j]->getFrameId()] = MapObservation(it_matches_cur_frame_j->second, local_map_frames[local_frame_j]);
          }
        }
        // Mark that feature has been checked
        b_feat_searched[feat_id] = true;

        // we need at least one match with feature in the focus frame
        if (feat_obs_nviews.size() < 1)
          continue;

        // Information for local frame
        Frame * frame_i;
        IntrinsicBase * cam_intrinsic_i;
        Mat3 K_inv_i;
        Mat3 frame_i_Rwc;

        // Test if matched feature correspondences for a new landmark
        // Get location of feature in current frame
        const Vec2 & pt_frame = frame->getCamCalibrated() ? frame->getFeaturePosition(feat_id) : cam_intrinsic->remove_disto(frame->getFeaturePosition(feat_id));

        Vec3 pt_3D;

        // Go through all matches of the point and find the one with biggest parallax
        MapObservation * mo_useful_ptr = nullptr;
        // And check if there is a pair with sufficient parallax
        for (LandmarkObservations::iterator iter_other_obs = feat_obs_nviews.begin(); iter_other_obs != feat_obs_nviews.end(); ++iter_other_obs)
        {
            MapObservation & m_o = iter_other_obs->second;
            // Information about the camera
            frame_i = m_o.frame_ptr;
            frame_i_Rwc = frame_i->getRotationMatrixInverse();
            if (cam_intrinsic_i != frame_i->getCameraIntrinsics())
            {
              cam_intrinsic_i = frame_i->getCameraIntrinsics();
              K_inv_i = frame_i->getK_inv();
            }

            // Get measurement in the i-th view
            const Vec2 & pt_frame_i = frame_i->getCamCalibrated() ? frame_i->getFeaturePosition(m_o.feat_id) : cam_intrinsic_i->remove_disto(frame_i->getFeaturePosition(m_o.feat_id));

            // Get projection matrix of second camera
            //TODO: How to use if we have relative cameras
            P_i = frame_i->getProjectionMatrix();

            // Triangulate results
            TriangulateDLT(P, pt_frame, P_i, pt_frame_i, &pt_3D);

            // Get projection on first image
            Vec3 pt_3D_frame;
            PoseEstimator::getRelativePointPosition(pt_3D,nullptr,pt_3D_frame,frame);

            // Project on second image
            Vec3 pt_3D_frame_i;
            PoseEstimator::getRelativePointPosition(pt_3D,nullptr,pt_3D_frame_i,frame_i);

            // Compute parallax angle between points
            //const Vec3 ray_cur = Vec3(frame_cur_Rwc * Vec3( K_inv_cur * Vec3( pt_cur( 0 ), pt_cur( 1 ), 1.0 ) )).normalized();
            //const Vec3 ray_i = Vec3(frame_i_Rwc * Vec3( K_inv_i * Vec3( pt_i( 0 ), pt_i( 1 ), 1.0 ) )).normalized();
            const Vec3 ray_frame = Vec3(frame_Rwc * pt_3D_frame).normalized();
            const Vec3 ray_frame_i = Vec3(frame_i_Rwc * pt_3D_frame_i).normalized();

            const double mag = ray_frame.norm() * ray_frame_i.norm();
            const double dotAngle = ray_frame.dot( ray_frame_i );
            const double cosParallax = clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 );
            // if parallax between two points is big enough we consider the point (cosParallax smaller than threshold)
            if (cosParallax < init_min_cos_angle_pt)
            {
              // Check if the points are inside the scene borders
              if (!frame->checkLandmarkPosition(pt_3D_frame))
              {
                break;
              }
              if (!frame_i->checkLandmarkPosition(pt_3D_frame_i))
              {
                break;
              }

              // Save the observation which can be used for initial triangulation
              mo_useful_ptr = & m_o;
              break;
            }
        }

        // If we didnt find an appropriate observation we skip point
        if (!mo_useful_ptr)
          continue;

        // Add observation from current frame
        feat_obs_nviews[frame->getFrameId()] = MapObservation(feat_id,frame);

        // Save number of all observations
        const size_t n_feat_obs_init = feat_obs_nviews.size();

        // Check in how many frames the triangulated points matches the observations
        Vec3 pt_3D_i;
        for (LandmarkObservations::iterator iter_other_obs = feat_obs_nviews.begin(); iter_other_obs != feat_obs_nviews.end();)
        {
          MapObservation & m_o = iter_other_obs->second;

          PoseEstimator::getRelativePointPosition(pt_3D,nullptr,pt_3D_i,m_o.frame_ptr);

          if (!m_o.frame_ptr->checkLandmarkPosition(pt_3D_i))
          {
            iter_other_obs = feat_obs_nviews.erase(iter_other_obs);
            continue;
          }

          if (!m_o.frame_ptr->checkFeatureAssociation(pt_3D_i,m_o.feat_id,5.991))
          {
            iter_other_obs = feat_obs_nviews.erase(iter_other_obs);
            continue;
          }
          iter_other_obs++;
        }

        // If we deleted any observations we still need at least 3 to add it
        if (feat_obs_nviews.size()!=n_feat_obs_init)
        {
          continue;
        }


        // Create a new landmark
        std::unique_ptr<MapLandmark> map_landmark_new = std::unique_ptr<MapLandmark>(new MapLandmark());
        map_landmark_new->X_ = pt_3D;
        // Set all observations
        map_landmark_new->obs_ = feat_obs_nviews;
        map_landmark_new->setNumberOfObservations(map_landmark_new->obs_.size());
        cartographer_->updateBestMapPointDescriptor(map_landmark_new.get());
        map_landmark_new->updateLastNormal(frame);
        map_landmark_new->updateNormal();

        // Mark this frame as the last frame that observed it
        map_landmark_new->last_local_map_frame_id_ = frame->getFrameId();
        // Mark as new triangulated point
        map_landmark_new->association_type_ = 6;

        // Add new triangulated point
        #pragma omp critical
        {
          vec_new_pts_3D.push_back(std::move(map_landmark_new));
        }

      }
    }
    vec_new_pts_3D.shrink_to_fit();

    std::cout<<"Tracker: [Track Triangulate] Concatenation and triangulation: "<<omp_get_wtime() - start_time<<" s\n";
  }


  double Tracker_Features::computeInlierRatioMatches(Frame * frame, Hash_Map<MapLandmark *,IndexT> & map_putative_matches_3D_pts_frame_idx)
  {
    if (map_putative_matches_3D_pts_frame_idx.size() == 0)
      return 1.0;
    // -------------------
    // -- Determine which of the points are actually inliers
    // -------------------
    size_t n_inliers = 0;
    const IndexT & frame_id = frame->getFrameId();
    IntrinsicBase * cam_intrinsic = frame->getCameraIntrinsics();
    const double & frame_reproj_thresh = frame->reproj_thresh_sq_;
    std::vector<MapLandmark*> & frame_map_points = frame->map_points_;

    for (auto iter_match: map_putative_matches_3D_pts_frame_idx)
    {
      // Get Landmark data
      const MapLandmark * map_landmark = iter_match.first;
      const IndexT feat_id_cur = iter_match.second;

      // Project point to frame coordinate system
      Vec3 pt_3D_frame;
      PoseEstimator::getRelativePointPosition(map_landmark->X_,map_landmark->ref_frame_,pt_3D_frame,frame);


      //if (frame->checkLandmarkPosition(pt_3D_frame) && frame->checkFeatureAssociation(pt_3D_frame,feat_id_cur,5.991))
      if (frame->checkFeatureAssociation(pt_3D_frame,feat_id_cur,5.991))
      {
        n_inliers++;
      }
    }
    return n_inliers / map_putative_matches_3D_pts_frame_idx.size();
  }


  void Tracker_Features::markInliersInFrame(Frame * frame, Hash_Map<MapLandmark *,IndexT> & map_putative_matches_3D_pts_frame_idx, bool b_check_error, size_t association_id)
  {
    // -------------------
    // -- Determine which of the points are actually inliers
    // -------------------
    const IndexT & frame_id = frame->getFrameId();
    IntrinsicBase * cam_intrinsic = frame->getCameraIntrinsics();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t pm_i = 0; pm_i < map_putative_matches_3D_pts_frame_idx.size(); ++pm_i)
    {
      Hash_Map<MapLandmark *,IndexT>::iterator iter_p_match = map_putative_matches_3D_pts_frame_idx.begin();
      std::advance(iter_p_match,pm_i);

      // Get Landmark data
      MapLandmark * map_landmark = iter_p_match->first;
      const IndexT feat_id_cur = iter_p_match->second;

      // Project point to frame coordinate system
      Vec3 pt_3D_frame;
      PoseEstimator::getRelativePointPosition(map_landmark->X_,map_landmark->ref_frame_,pt_3D_frame,frame);

      //if (!b_check_error || (frame->checkLandmarkPosition(pt_3D_frame) && frame->checkFeatureAssociation(pt_3D_frame,feat_id_cur,5.991)))
      if (!b_check_error || frame->checkFeatureAssociation(pt_3D_frame,feat_id_cur,5.991))
      {
        frame->setMapPoint(feat_id_cur,map_landmark);
        map_landmark->increaseNumberOfObservations();
        map_landmark->association_type_ = association_id;
      }
    }
  }

  void Tracker_Features::removeOutliersInFrame(Frame * frame)
  {
    // -------------------
    // -- Determine which of the points are actually inliers
    // -------------------
    const IndexT & frame_id = frame->getFrameId();
    IntrinsicBase * cam_intrinsic = frame->getCameraIntrinsics();
    const double & frame_reproj_thresh = frame->reproj_thresh_sq_;
    std::vector<MapLandmark*> & frame_map_points = frame->map_points_;

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t pm_i = 0; pm_i < frame_map_points.size(); ++pm_i)
    {
      bool b_inlier;
      if (!frame_map_points[pm_i])
        continue;

      MapLandmark * map_landmark = frame_map_points[pm_i];

      // Project point to frame coordinate system
      Vec3 pt_3D_frame;
      PoseEstimator::getRelativePointPosition(map_landmark->X_,map_landmark->ref_frame_,pt_3D_frame,frame);

      if (pt_3D_frame(2) < 0 || !frame->checkFeatureAssociation(pt_3D_frame,pm_i,5.991))
      //if (!frame->checkLandmarkPosition(pt_3D_frame) || !frame->checkFeatureAssociation(pt_3D_frame,pm_i,5.991))
      {
        frame->clearMapPoint(pm_i);
        map_landmark->decreaseNumberOfObservations();
      }
    }
  }



  void Tracker_Features::removeOutliersInNewTriangulatedPoints(Frame * frame, std::vector<std::unique_ptr<MapLandmark> > & vec_putative_new_pts_3D)
  {
    // Check newly triangulated points if they have at least two inliers
    std::vector<size_t> vec_new_pts_outliers;
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for(size_t t_i = 0; t_i < vec_putative_new_pts_3D.size(); ++t_i)
    {
      MapLandmark * map_landmark = vec_putative_new_pts_3D[t_i].get();

      LandmarkObservations & obs = map_landmark->obs_;
      const size_t n_feat_obs_init = obs.size();

      bool b_sufficient_parallax = false;

      // Get observation in current frame
      Vec3 pt_3D_frame;
      size_t frame_id = frame->getFrameId();
      const MapObservation & m_o = obs.find(frame_id)->second;
      // Check if everything is ok in current frame
      PoseEstimator::getRelativePointPosition(map_landmark->X_,map_landmark->ref_frame_,pt_3D_frame,frame);
      if (!frame->checkLandmarkPosition(pt_3D_frame) || !frame->checkFeatureAssociation(pt_3D_frame,m_o.feat_id,5.991))
      {
        continue;
      }

      // Check the rest of the frames and check if there is any point with sufficient parallax
      for(LandmarkObservations::iterator iter_mo = obs.begin(); iter_mo != obs.end();)
      {
        MapObservation & m_o_i =  iter_mo->second;
        Frame * & frame_i = m_o_i.frame_ptr;
        // Put 3D point into coordinate system of frame
        Vec3 pt_3D_frame_i;
        PoseEstimator::getRelativePointPosition(map_landmark->X_,map_landmark->ref_frame_,pt_3D_frame_i,frame_i);
        if (!frame_i->checkLandmarkPosition(pt_3D_frame_i) || !frame_i->checkFeatureAssociation(pt_3D_frame_i,m_o_i.feat_id,5.991))
        {
          // Remove measurement from 3D point
          iter_mo = obs.erase(iter_mo);
          map_landmark->decreaseNumberOfObservations();
        }
        else
        {
          if (!b_sufficient_parallax)
          {
            // Compute parallax angle
            const Vec3 ray_frame = Vec3(frame->getRotationMatrixInverse() * pt_3D_frame).normalized();
            const Vec3 ray_frame_i = Vec3(frame_i->getRotationMatrixInverse() * pt_3D_frame_i).normalized();

            const double mag = ray_frame.norm() * ray_frame_i.norm();
            const double dotAngle = ray_frame.dot( ray_frame_i );
            const double cosParallax = clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 );

            if (cosParallax < init_min_cos_angle_pt)
            {
              b_sufficient_parallax = true;
            }
          }
          ++iter_mo;
        }
      }

      if ((obs.size()!=n_feat_obs_init) || !b_sufficient_parallax)
      {
        // Delete point from the vector of new points
        #pragma omp critical
        {
          vec_new_pts_outliers.push_back(t_i);
        }
      }
      else if (obs.size()>2)
      {
        // Mark that point was triangulated from multiple points
        map_landmark->association_type_ = 7;
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
  void Tracker_Features::detect
  (
    const image::Image<unsigned char> & ima,
    Frame * frame,
    const size_t & min_count,
    const image::Image<unsigned char> * mask
  )
  {
    // Detect feature points
    size_t n_feats_detected = featureExtractor_->detect(ima,frame,min_count,mask);

    if (!(n_feats_detected > 0))
    {
      return;
    }

    // Describe detected features
    featureExtractor_->describe(ima,frame);
    // Undistort points
    frame->updateFeaturesData();
  }

} // namespace VO
} // namespace openMVG
