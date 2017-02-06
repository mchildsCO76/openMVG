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
#include <openMVG/vsslam/optimization/Optimization.hpp>
#include "openMVG/multiview/triangulation.hpp"

namespace openMVG  {
namespace VSSLAM  {

class Tracker_Features : public Abstract_Tracker
{
public:
  // ---------------
  // Parameters
  // ---------------

  /// Feature extractor
  Abstract_FeatureExtractor * featureExtractor_;

  // Tracking data
  size_t max_tracked_points = 1500;
  Hash_Map<size_t,MapLandmark *> tracking_cur_idx_3D_ptr;

  Hash_Map<size_t,std::vector<MapObservation *> > new_feat_cur_ref_ids;

  // Tracking settings
  size_t min_init_ref_tracks = 500; // Min number of tracks reference frame for initialization has to have
  size_t min_frame_tracks = 300; // Min number of tracks frame has to have
  size_t max_frame_tracks = 0; // Max number of feats detected in a frame (0 - unlimited)
  size_t min_matches_init_pose = 500; // Min number of matches for init pose estimation


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

  void clearTrackingData()
  {
    new_feat_cur_ref_ids.clear();
    tracking_cur_idx_3D_ptr.clear();
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
    // Clear data for tracking from prev frame
    clearTrackingData();

    // Detect features
    std::cout<<"Start detection:\n";
    double startTime = omp_get_wtime();

    detect(ima,mCurrentFrame,max_tracked_points,0);

    size_t n_feats_detected = mCurrentFrame->N_feats();
    std::cout<<"Candidate features: "<<n_feats_detected<<"\n";

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Detect features time:"<<secsElapsed<<"\n";


    // If tracking is not initialized
    if (trackingStatus == TRACKING_STATUS::NOT_INIT)
    {
      // Detect new features and add frame as init reference frame
      std::cout<<"DETECT FROM STRACH A!\n";

      // Check if enough features are detected
      if (n_feats_detected > min_init_ref_tracks)
      {
        track_status = trySystemInitialization(ima);
      }
      else
      {
        // Insuccessful detection of features
        std::cout<<"Insufficient number of features detected!\n";
        track_status = false;
      }
    }
    else if (trackingStatus == TRACKING_STATUS::INIT)
    {
      std::cout<<"TRY TO TRACK FROM INIT REF FRAME ID: "<<init_ref_frame->frameId_<<": features N: "<<init_ref_frame->N_feats()<<"!\n";

      // Check if enough features are detected
      if (n_feats_detected > min_init_ref_tracks)
      {
        // Try initializing the system
        track_status = trySystemInitialization(ima);
      }
      else
      {
        // not enough features detected
        resetSystemInitialization();
        trackingStatus = Abstract_Tracker::TRACKING_STATUS::NOT_INIT;
        track_status = false;
      }
    }
    else if (trackingStatus == TRACKING_STATUS::OK)
    {
      size_t min_track_matches = 30;
      bool bSuccessfulTrack = false;

      // ----------------
      // -- Motion Model tracking
      // ----------------
      // Use previous tracking to predict where the pose will be and try matching with previously triangulated map points
      if (motionModel.isValid())
      {
        std::cout<<"Motion model is valid\n";

        // Predict location of current frame
        Mat4 T_predict = motionModel.predictLocation(mPrevFrame.get(),mCurrentFrame.get());
        mCurrentFrame->setPose_T(T_predict);

        // Match current frame with lastFrame
        bSuccessfulTrack = trackWithMotionModel(min_track_matches);
      }
      std::cout<<"TRACK MM: "<<bSuccessfulTrack<<"\n";

      // ----------------
      // -- Reference frame tracking
      // ----------------
      // If motion model tracking didnt work we try with reference frame
      if (!bSuccessfulTrack)
      {
        // Try with normal keyframe tracking
        bSuccessfulTrack = trackWithReferenceFrame(min_track_matches);
      }

      // ----------------
      // -- Local map tracking
      // ----------------
      if (bSuccessfulTrack)
      {
        // Try to track other points from local map

        trackLocalMap();

        // Do pose BA

        // Decide if its a keyframe
        bool bKeyframe = true;

        // Add current frame as new keyframe
        if (bKeyframe)
        {
          cartographer_->addKeyFrameToMap(mCurrentFrame);
          cartographer_->addObservations(mCurrentFrame.get());
          // Set current frame as the last keyframe
          mCurrentFrame->updateFrameVisibilityConnections();
          mLastRefFrame = mCurrentFrame->share_ptr();
        }

        // Set motion model
        motionModel.updateMotionModel(mPrevFrame.get(),mCurrentFrame.get());

        track_status = true;
      }
      else
      {
        // We couldnt match with previous frame
        if (cartographer_->getNumberOfKeyFrames() < 5)
        {
          std::cout<<"TRACKING: System lost before 5 frames! :: New init required\n";
          // Set tracking as not initialized
          resetSystemInitialization();
          trackingStatus = Abstract_Tracker::TRACKING_STATUS::NOT_INIT;
          track_status = false;
        }
        else
        {
          trackingStatus = Abstract_Tracker::TRACKING_STATUS::LOST;
          track_status = false;
        }
      }

    }

    std::cout<<"Tracks available: "<<mCurrentFrame->N_feats()<<"\n";

    mPrevFrame.swap(mCurrentFrame);
    mCurrentFrame.reset();

    // Return if tracking is ok
    return track_status;
  }


  /// INITIALIZATION

  void resetSystemInitialization()
  {
    std::cout<<"Reset system initialization!\n";
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
    std::cout<<"Set new reference initialization frame\n";
  }

  bool trySystemInitialization
  (
    const image::Image<unsigned char> & ima
  )
  {
    std::cout<<"System initialized process B!\n";

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
        setReferenceSystemInitialization(mCurrentFrame);
        // Set system to INIT
        trackingStatus = Abstract_Tracker::TRACKING_STATUS::INIT;
        return true;
      }

      // -------------------
      // -- Match current frame to init reference frame (all-all 2D-2D)
      // -------------------
      // Matching settings
      float desc_ratio = 0.8;

      Hash_Map<size_t,size_t> putative_matches_ref_cur_idx;

      double startTime = omp_get_wtime();

      size_t n_put_matches = matching_AllAll_2D_2D(featureExtractor_, init_ref_frame.get(), mCurrentFrame.get(), putative_matches_ref_cur_idx, desc_ratio);

      double stopTime = omp_get_wtime();
      double secsElapsed = stopTime - startTime; // that's all !
      std::cout<<"Matching time (no MM):"<<secsElapsed<<"\n";

      std::cout<<"Matches with ref frame: "<<n_put_matches<<"\n";


      // Check if enough matches with reference image
      if (n_put_matches < min_matches_init_pose)
      {
        // Not enough matches with initialization reference frame
        // Set this frame as initialization reference frame and try with next frame
        std::cout<<"Not enough matches for initialization process - setting this as REF frame\n";
        setReferenceSystemInitialization(mCurrentFrame);
        trackingStatus = Abstract_Tracker::TRACKING_STATUS::INIT;
        return true;
      }

      // -------------------
      // -- Find relative pose cur_T_ref = [ cur_R_ref t_ref]
      // -------------------

      // Info about the frames
      const size_t frame_id_ref = init_ref_frame->frameId_;
      const size_t frame_id_cur = mCurrentFrame->frameId_;
      const Camera * cam_ref = init_ref_frame->cam_;
      const Camera * cam_cur = mCurrentFrame->cam_;
      const IntrinsicBase * cam_intrinsic_ref = cam_ref->cam_intrinsic_ptr;
      const IntrinsicBase * cam_intrinsic_cur = cam_cur->cam_intrinsic_ptr;


      bool bSuccessMotion = false;
      // Min angle between rays for the point to be triangulated (0.99995 ~ 0.5deg; 0.99998 ~ 0.36deg; 0.9998 ~ 1.15deg)
      double min_cos_angle_init_point = 0.99995;
      // Min points needed for successful initialization
      size_t min_init_map_points = 50;

      // Recovered motion
      Mat3 c_R_r;
      Vec3 c_t_r; // Center of origin of w (reference cam) in c (current cam)
      double c_s_r;
      double AC_reprojection_thresh;
      std::vector<size_t> inliers_M;  // Inliers of an estimated model

      startTime = omp_get_wtime();

      // If calibrated camera we use Essential matrix otherwise Fundamental
      if (cam_ref->bCalibrated && cam_cur->bCalibrated)
      {
        bSuccessMotion = estimateRobustRelativePoseHE(init_ref_frame, mCurrentFrame, putative_matches_ref_cur_idx, min_cos_angle_init_point, min_init_map_points, c_R_r, c_t_r, inliers_M, AC_reprojection_thresh);
        // Estimate scale of current camera
        c_s_r = 1.0;
      }
      else
      {
        // TODO: Compute HF
      }

      stopTime = omp_get_wtime();
      secsElapsed = stopTime - startTime; // that's all !
      std::cout<<"Motion  time:"<<secsElapsed<<"\n";

      if(!bSuccessMotion)
      {
        std::cout<<"TRACKING: Motion estimation failed\n";
        return false;
      }
      // -------------------
      // -- Set initial pose estimations
      // -------------------
      init_ref_frame->setReferenceFrame_Rts(nullptr,Mat3::Identity(),Vec3::Zero(),1.0);
      init_ref_frame->AC_reprojection_thresh_ = AC_reprojection_thresh;
      std::cout<<"F1 T: "<<init_ref_frame->AC_reprojection_thresh_ <<"\n";

      mCurrentFrame->setReferenceFrame_Rts(nullptr,c_R_r,c_t_r,c_s_r);
      mCurrentFrame->AC_reprojection_thresh_ = AC_reprojection_thresh;
      std::cout<<"F2 T: "<<mCurrentFrame->AC_reprojection_thresh_ <<"\n";


      // -------------------
      // -- Triangulate inliers
      // -------------------
      std::vector<Frame*> local_frames(2);
      std::vector<Vec3> triangulated_pts(inliers_M.size());
      std::vector<std::pair<Vec3, std::deque<std::pair<Frame*,size_t> > > > vec_new_triangulated_pts(inliers_M.size());

      // Add local frames
      local_frames[0] = init_ref_frame.get();
      local_frames[1] = mCurrentFrame.get();

      // Estimate projection matrix of both cameras based on estimated poses
      Mat34 P1,P2;
      switch (cartographer_->getMapPointType())
      {
        case VSSLAM::MAP_CAMERA_TYPE::ABSOLUTE:
          P1 = init_ref_frame->getProjectionMatrix();
          P2 = mCurrentFrame->getProjectionMatrix();
        break;
        case VSSLAM::MAP_CAMERA_TYPE::RELATIVE:
        /*{
          // If cameras are relative we find tranformation
          Mat4 T_cam1_2;
          if (!getRelativeCameraTransformation(frame_1.get(), frame_2.get(),T_cam1_2))
          {
            std::cout<<"Cant find relative poses!!\n Something is seriously wrong!\n";
            return false;
          }
          P1 = frame_1->getK()*T_cam1_2.block(0,0,3,4);
          // Set second camera as base
          P2 = HStack(Mat3::Identity(), Vec3::Zero());
        }*/
        break;
      }

      // Iterate through inliers and triangulate
      for ( size_t k_i = 0; k_i < inliers_M.size(); ++k_i)
      {
        Hash_Map<size_t,size_t>::iterator m_iter = putative_matches_ref_cur_idx.begin();
        std::advance(m_iter,inliers_M[k_i]);

        std::deque<std::pair<Frame*,size_t> > & measurements = vec_new_triangulated_pts[k_i].second;

        const size_t feat_id_ref = m_iter->first;
        const size_t feat_id_cur = m_iter->second;

        // Position of detected features
        //   If we calibrated camera we just use undistorted
        //   If not calibrated we remove the distortion as best as we know it (current estimation)
        const Vec2
          & x1_ = cam_ref->bCalibrated ? init_ref_frame->getFeaturePositionUndistorted(feat_id_ref) : cam_intrinsic_ref->remove_disto(init_ref_frame->getFeaturePositionUndistorted(feat_id_ref)),
          & x2_ = cam_cur->bCalibrated ? mCurrentFrame->getFeaturePositionUndistorted(feat_id_cur) : cam_intrinsic_cur->remove_disto(mCurrentFrame->getFeaturePositionUndistorted(feat_id_cur));


        // Triangulate results
        TriangulateDLT(P1, x1_, P2, x2_, &(vec_new_triangulated_pts[k_i].first));

        // Add observations
        measurements.push_back(std::make_pair(init_ref_frame.get(),feat_id_ref));
        measurements.push_back(std::make_pair(mCurrentFrame.get(),feat_id_cur));
      }

      // -------------------
      // -- Optimize initial poses
      // -------------------
      Hash_Map<MapLandmark *,size_t> matches_3D_ptr_cur_idx;
      bool bSuccessInitPoseOptimization = false;

      VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
      options.linear_solver_type_ = ceres::DENSE_SCHUR;
      VSSLAM_Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
      if (!bundle_adjustment_obj.OptimizePose(local_frames, nullptr, &vec_new_triangulated_pts))
      {
        return false;
      }

      // -------------------
      // -- Determine inliers
      // -------------------
      std::vector<std::pair<Vec3, std::deque<std::pair<Frame*,size_t> > > > vec_new_pts_3D_obs;

      Similarity3 & pose_ref = init_ref_frame->pose_;
      Mat3 cam_ref_Rcw = init_ref_frame->getInverseRotationMatrix();

      Similarity3 & pose_cur = mCurrentFrame->pose_;
      Mat3 cam_cur_Rcw = mCurrentFrame->getInverseRotationMatrix();

      // Go through all the points and check if they are inliers
      for(auto point_it: vec_new_triangulated_pts)
      {
        Vec3 & pt3D = point_it.first;
        std::deque<std::pair<Frame*,size_t> > & measurements = point_it.second;

        // Check points for angle and reprojection error
        const size_t feat_ref_id = measurements[0].second;
        const size_t feat_cur_id = measurements[1].second;

        // Check that the point is infront of the camera
        if (pose_ref.depth(pt3D) <= 0 || pose_cur.depth(pt3D) <= 0)
          continue;

        const Vec2 & pt_ref = init_ref_frame->getFeaturePositionUndistorted(feat_ref_id);
        const Vec2 & pt_cur = mCurrentFrame->getFeaturePositionUndistorted(feat_cur_id);

        // Ray between point and focal center
        const Vec3 ray_ref = Vec3(cam_ref_Rcw * cam_intrinsic_ref->operator ()(pt_ref)).normalized();
        const Vec3 ray_cur = Vec3(cam_cur_Rcw * cam_intrinsic_cur->operator ()(pt_cur)).normalized();

        const double mag = ray_ref.norm() * ray_cur.norm();
        const double dotAngle = ray_ref.dot( ray_cur );
        const double cosParallax = clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 );

        if (cosParallax > min_cos_angle_init_point)
        {
          continue;
        }

        const Vec2 pt_proj_ref = cam_intrinsic_ref->cam2ima(cam_intrinsic_ref->have_disto()?cam_intrinsic_ref->add_disto(Vec3(pose_ref(pt3D)).hnormalized()):Vec3(pose_ref(pt3D)).hnormalized());
        if ((pt_ref - pt_proj_ref).squaredNorm() > init_ref_frame->AC_reprojection_thresh_ )
          continue;

        const Vec2 pt_proj_cur = cam_intrinsic_cur->cam2ima(cam_intrinsic_cur->have_disto()?cam_intrinsic_cur->add_disto(Vec3(pose_cur(pt3D)).hnormalized()):Vec3(pose_cur(pt3D)).hnormalized());
        if ((pt_cur - pt_proj_cur).squaredNorm() > mCurrentFrame->AC_reprojection_thresh_ )
          continue;

        vec_new_pts_3D_obs.push_back(std::make_pair(pt3D,measurements));
      }

        /*
        // Computed rays
        std::deque<std::pair<Vec3,std::pair<Frame*,size_t>*> > pt_rays;

        for(auto & measurement : point_it.second)
        {
          Frame * frame_i = measurement.first;
          const IntrinsicBase * cam_intrinsic_i = frame_i->cam_->cam_intrinsic_ptr;
          Similarity3 & pose_i = frame_i->pose_;
          Mat3 cam_i_Rcw = frame_i->getInverseRotationMatrix();

          const size_t feat_i_id = measurement.second;

          // Check that the point is infront of the camera
          if (pose_i.depth(pt3D) <= 0)
            continue;

          const Vec2 & pt = frame_i->getFeaturePositionUndistorted(feat_i_id);

          // Compute residual error in frame 1
          // const Vec2 pt_1 = cam_intrinsic_1->have_disto()? cam_intrinsic_1->cam2ima(cam_intrinsic_1->add_disto(Vec3(pose_1(landmark.X_)).hnormalized())) :  cam_intrinsic_1->cam2ima(Vec3(pose_1(landmark.X_)).hnormalized());
          const Vec2 pt_proj = cam_intrinsic_i->cam2ima(cam_intrinsic_i->have_disto()?cam_intrinsic_i->add_disto(Vec3(pose_i(pt3D)).hnormalized()):Vec3(pose_i(pt3D)).hnormalized());
          //std::cout<<"AA: "<<frame_i->AC_reprojection_thresh_<<" BB: "<<(pt - pt_proj).squaredNorm()<<" CC: "<<pt<<" DD: "<<pt_proj<<" EE: "<<pt3D<<"\n";
          if ((pt - pt_proj).squaredNorm() > frame_i->AC_reprojection_thresh_)
            continue;

          // Ray between point and focal center
          const Vec3 ray1 = Vec3(cam_i_Rcw * cam_intrinsic_i->operator ()(pt)).normalized();
          pt_rays.emplace_back(std::make_pair(ray1,&measurement));
        }

        // Measurements that have at least one point with sufficient angle
        std::deque<std::pair<Frame*,size_t> > accepted_measurements;

        // Check which ray has at least one that is far enough
        for (size_t r_i = 0; r_i < pt_rays.size(); ++r_i)
        {
          for (size_t r_j = 0; r_j < pt_rays.size(); ++r_j)
          {
            if (r_i != r_j)
            {
              const double mag = pt_rays[r_i].first.norm() * pt_rays[r_j].first.norm();
              const double dotAngle = pt_rays[r_i].first.dot( pt_rays[r_j].first );
              const double cosParallax = clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 );

              if (cosParallax < min_cos_angle_init_point)
              {
                accepted_measurements.push_back(*(pt_rays[r_i].second));
                break;
              }
            }
          }
        }
        // At least two measurements for triangulation
        if (accepted_measurements.size() > 1)
        {
          vec_new_pts_3D_obs.push_back(std::make_pair(pt3D,accepted_measurements));
        }
      }*/

      // -------------------
      // -- Initialize Map
      // -------------------
      cartographer_->initializeMap(init_ref_frame,mCurrentFrame, vec_new_pts_3D_obs);


      // Save inital map
      Save_PLY(cartographer_->slam_data,"Initial_Map.ply",sfm::ESfM_Data::ALL);

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

  void triangulateNewPoints()
  {

  }


  void trackLocalMap()
  {
    std::cout<<"Track local map\n";
    size_t win_size = (20*20);
    float desc_ratio = 0.8;
    // Generate local map
    std::vector<Frame *> local_map_frames;
    std::vector<MapLandmark *> local_map_points;

    // Get frames that see points that current frame sees
    mCurrentFrame->getFrameVisibilityConnections(local_map_frames);
    std::cout<<"Frame visibility: "<<local_map_frames.size()<<"\n";
    // Get all unmatched map points from neighbor frames
    local_map_points = cartographer_->getLocalMapPoints(mCurrentFrame.get(),local_map_frames);
    std::cout<<"points visibility: "<<local_map_points.size()<<"\n";
    std::cout<<"Local map: "<<local_map_frames.size()<<" :: "<<local_map_points.size()<<"\n";

    // Try matching with the rest of the frames
    //matchingByProjection(featureExtractor_, local_map_points, mCurrentFrame.get(), matches_3D_ptr_cur_idx, win_size, desc_ratio);

    //std::cout<<"Matces with local map: "<<matches_3D_ptr_cur_idx.size()<<"\n";

  }


  bool trackWithMotionModel(const size_t & min_track_matches)
  {
    // Matching settings
    size_t win_size = (20*20);
    float desc_ratio = 0.8;

    Hash_Map<MapLandmark *,size_t> putative_matches_3D_ptr_cur_idx;

    // Try matching with motion model
    matchingByProjection(featureExtractor_, mPrevFrame->map_points_, mCurrentFrame.get(), putative_matches_3D_ptr_cur_idx, win_size, desc_ratio);

    std::cout<<"Matched by projecton: "<<putative_matches_3D_ptr_cur_idx.size()<<"\n";

    if (putative_matches_3D_ptr_cur_idx.size() < min_track_matches)
    {
      std::cout<<"Matching with MM failed - try 1!\n";
      // Try matching with projection with larger window
      matchingByProjection(featureExtractor_, mPrevFrame->map_points_, mCurrentFrame.get(), putative_matches_3D_ptr_cur_idx, win_size*2, desc_ratio);
      if (putative_matches_3D_ptr_cur_idx.size() < min_track_matches)
      {
        std::cout<<"Matching with MM failed - try 2!\n";
        return false;
      }
    }

    // Show how many matches we actually found
    std::cout<<"Matches by projecton OK: "<<putative_matches_3D_ptr_cur_idx.size()<<"\n";


    // -------------------
    // -- Do tiny BA to refine the pose
    // --- refine only Rotations & translations & scale (keep intrinsic constant)
    // -------------------
    VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
    VSSLAM_Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    std::vector<Frame*> tmp = {mCurrentFrame.get()};
    if (!bundle_adjustment_obj.OptimizePose(tmp, &putative_matches_3D_ptr_cur_idx, nullptr))
    {
      return false;
    }
    /*
    if (!bundle_adjustment_obj.OptimizePose_2D_3D(mCurrentFrame.get(),putative_matches_3D_ptr_cur_idx))
    {
      return false;
    }*/

    // TODO: Check how to set reprojection threshold if we just use BA
    mCurrentFrame->AC_reprojection_thresh_ = 16.0;

    // Determine which of the points are inliers
    const Similarity3 & pose_cur = mCurrentFrame->pose_;
    const IntrinsicBase * cam_cur_intrinsic = mCurrentFrame->cam_->cam_intrinsic_ptr;
    const double & frame_reproj_error = mCurrentFrame->AC_reprojection_thresh_;

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t pm_i = 0; pm_i < putative_matches_3D_ptr_cur_idx.size(); ++pm_i)
    {
      Hash_Map<MapLandmark *,size_t>::iterator iter_p_match = putative_matches_3D_ptr_cur_idx.begin();
      std::advance(iter_p_match,pm_i);
      // Get Landmark data
      MapLandmark * map_point = iter_p_match->first;
      size_t feat_id_cur = iter_p_match->second;

      // Check that the point is infront of the camera
      if (pose_cur.depth(map_point->X_) <= 0)
      {
        continue;
      }

      // Get the observation
      const Vec2 & obs_cur = mCurrentFrame->getFeaturePositionUndistorted(feat_id_cur);

      // Compute residual error in current frame
      // We add distortion to 3D points if we have it
      const Vec2 pt_cur = cam_cur_intrinsic->cam2ima(cam_cur_intrinsic->have_disto()?cam_cur_intrinsic->add_disto(Vec3(pose_cur(map_point->X_)).hnormalized()):Vec3(pose_cur(map_point->X_)).hnormalized());

      if ((obs_cur - pt_cur).squaredNorm() < frame_reproj_error)
      {
        mCurrentFrame->map_points_[feat_id_cur] = map_point;
      }
    }

    size_t n_matches = 0;
    for (size_t i = 0; i<mCurrentFrame->map_points_.size(); ++i)
    {
      if (mCurrentFrame->map_points_[i])
        n_matches++;
    }

    std::cout<<"Matches after BA: "<<n_matches<<"\n";


    if (n_matches < min_track_matches)
      return false;
    return true;
  }



  bool trackWithReferenceFrame(const size_t & min_track_matches)
  {
    // Matching settings
    float desc_ratio = 0.8;

    Hash_Map<MapLandmark* ,size_t> putative_matches_3D_ptr_cur_idx;

    matchingAllAll_3D_2D(featureExtractor_, mLastRefFrame->map_points_, mCurrentFrame.get(), putative_matches_3D_ptr_cur_idx, desc_ratio);

    std::cout<<"Matched AllAll with "<<mLastRefFrame->frameId_<<" last ref: "<<putative_matches_3D_ptr_cur_idx.size()<<"\n";

    // Try matching with motion model
    if (putative_matches_3D_ptr_cur_idx.size() < min_track_matches)
    {
      std::cout<<"Matching AllAll failed! Only :"<<putative_matches_3D_ptr_cur_idx.size()<<" matches!\n";
      return false;
    }

    // Estimate the pose using AC P3P
    // Perform resection to get the camera pose
    Mat3 R;
    Vec3 t;
    double s;
    std::vector<size_t> vec_inliers;
    double AC_threshold;
    // TODO: How to estimate scale!!
    bool bSuccessMotion = estimateRobustPose(mCurrentFrame.get(),putative_matches_3D_ptr_cur_idx, min_track_matches,R,t,vec_inliers,AC_threshold);
    s = 1.0f;

    if (!bSuccessMotion)
    {
      std::cout<<"Robust pose estimation with AllAll failed!\n";
      return false;
    }

    // Update estimated pose
    mCurrentFrame->setReferenceFrame_Rts(nullptr,R,t,s);
    mCurrentFrame->AC_reprojection_thresh_ = AC_threshold;

    // Mark inliers in the set
    //#ifdef OPENMVG_USE_OPENMP
    //#pragma omp parallel for schedule(dynamic)
    //#endif
    for (size_t k: vec_inliers)
    {
      Hash_Map<MapLandmark *, size_t>::iterator it_p_match = putative_matches_3D_ptr_cur_idx.begin();
      std::advance(it_p_match,k);

      mCurrentFrame->map_points_[it_p_match->second] = it_p_match->first;
    }


    size_t n_matches = 0;
    for (size_t i = 0; i<mCurrentFrame->map_points_.size(); ++i)
    {
      if (mCurrentFrame->map_points_[i])
        n_matches++;
    }

    std::cout<<"Matches all-all: "<<n_matches<<"\n";

    if (n_matches < min_track_matches)
      return false;
    return true;
  }





  bool detect
  (
    const image::Image<unsigned char> & ima,
    std::shared_ptr<Frame> & frame,
    const size_t min_count,
    const size_t max_count
  )
  {
    // Detect feature points
    //double startTime = omp_get_wtime();

    size_t n_feats_detected = featureExtractor_->detect(ima,frame->regions_,min_count,max_count);

    //double stopTime = omp_get_wtime();
    //double secsElapsed = stopTime - startTime; // that's all !
    //std::cout<<"Detect time:"<<secsElapsed<<"\n";
    // Describe detected features

    //startTime = omp_get_wtime();

    featureExtractor_->describe(ima,frame->regions_.get());

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
