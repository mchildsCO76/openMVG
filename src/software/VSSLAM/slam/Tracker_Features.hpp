// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>

#include <openMVG/types.hpp>
#include <software/VSSLAM/slam/Frame.hpp>
#include <software/VSSLAM/slam/Abstract_Tracker.hpp>
#include <software/VSSLAM/slam/Abstract_FeatureExtractor.hpp>

namespace openMVG  {
namespace VSSLAM  {

struct Tracker_Features : public Abstract_Tracker
{
  /// Feature extractor
  Abstract_FeatureExtractor * featureExtractor_;

  // Tracking candidates
  //std::unique_ptr<features::Regions> candidate_regions;
  //features::PointFeatures candidate_pts_undist;
  //std::vector<bool> candidate_pts_used;

  // Tracking data
  Hash_Map<size_t,size_t> tracking_feat_cur_ref_ids;

  // Tracking settings
  size_t min_init_ref_tracks = 500; // Min number of tracks reference frame for initialization has to have
  size_t min_frame_tracks = 300; // Min number of tracks frame has to have
  size_t max_frame_tracks = 0; // Max number of feats detected in a frame (0 - unlimited)
  size_t min_matches_init_pose = 500; // Min number of matches for init pose estimation



  Tracker_Features(){}

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
    tracking_feat_cur_ref_ids.clear();

  }

  /// Try to track current point set in the provided image
  /// return false when tracking failed (=> to send frame to relocalization)
  bool track
  (
    const image::Image<unsigned char> & ima,
    std::shared_ptr<Frame> current_frame
  ) override
  {
    // Set current frame
    mCurrentFrame = current_frame->share_ptr();
    // Clear data for tracking from prev frame
    tracking_feat_cur_ref_ids.clear();


    // Detect features
    std::cout<<"Start detection:\n";
    double startTime = omp_get_wtime();

    detect(ima,mCurrentFrame,max_tracked_points,0);
    size_t n_feats_detected = mCurrentFrame->regions->RegionCount();
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
        trySystemInitialization(ima);
      }
      else
      {
        // Insuccessful detection of features
        std::cout<<"Insufficient number of features detected!\n";
      }
    }
    else if (trackingStatus == TRACKING_STATUS::INIT)
    {
      std::cout<<"TRY TO TRACK FROM INIT REF FRAME ID: "<<init_ref_frame->frameId_<<": features"<<init_ref_frame->getTracksSize()<<"!\n";

      // Check if enough features are detected
      if (n_feats_detected > min_init_ref_tracks)
      {
        // Try initializing the system
        trySystemInitialization(ima);
      }
      else
      {
        // not enough features detected
        resetSystemInitialization();
      }

    }


    std::cout<<"Tracks available: "<<mCurrentFrame->regions->RegionCount()<<"\n";

    mPrevFrame.swap(mCurrentFrame);
    mCurrentFrame.reset();

    // Return if tracking is ok
    return trackingStatus==TRACKING_STATUS::OK;
  }
/*
  // suggest new feature point for tracking (count point are kept)
  bool detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count,
    const size_t min_count
  ) const override
  {
    size_t n_feats = featureExtractor_->detect(ima,pt_to_track,count);

    if (n_feats < min_count)
      return true;
    return true;
  }
*/
  /// INITIALIZATION

  void resetSystemInitialization()
  {
    std::cout<<"Reset system initialization!\n";
    if (init_ref_frame)
    {
      init_ref_frame.reset();
    }
    trackingStatus = Abstract_Tracker::TRACKING_STATUS::NOT_INIT;
  }

  void setReferenceSystemInitialization(std::shared_ptr<Frame> & frame)
  {
    std::cout<<"System initialized process A!\n";

    // Reset all initialization settings
    resetSystemInitialization();
    // Set current frame as the new reference frame for initialization
    init_ref_frame = frame->share_ptr();

    trackingStatus = Abstract_Tracker::TRACKING_STATUS::INIT;
    std::cout<<"Set new reference initialization frame\n";
  }

  void trySystemInitialization
  (
    const image::Image<unsigned char> & ima
  )
  {
    std::cout<<"System initialized process B!\n";

    if (trackingStatus == Abstract_Tracker::TRACKING_STATUS::NOT_INIT)
    {
      // Clear all data that could be initialized

      // Set new reference initialization frame
      init_ref_frame = mCurrentFrame->share_ptr();

      // Set system to INIT
      trackingStatus = Abstract_Tracker::TRACKING_STATUS::INIT;
    }
    else if (trackingStatus == Abstract_Tracker::TRACKING_STATUS::INIT)
    {
      // Check that we have actual reference image
      if (!init_ref_frame)
      {
        resetSystemInitialization();
        return;
      }

      // -------------------
      // -- Match current frame to init reference frame (no MM)
      // -------------------

      // Matching settings
      size_t win_size = 50;
      float desc_ratio = 0.8;

      // -------------------
      // -- Find matches reference init frame - currrent frame
      // -------------------
      std::cout<<"Try matching with init reference image\n";
      Hash_Map<size_t,size_t> match_cur_ref_idx;

      double startTime = omp_get_wtime();

      size_t n_feat_frame_matches = matchFramesFeatureMatchingNoMM(init_ref_frame->regions.get(), init_ref_frame->pts_undist, mCurrentFrame->regions.get(), mCurrentFrame->pts_undist, match_cur_ref_idx, win_size, desc_ratio);

      std::cout<<"Matches with ref frame: "<<n_feat_frame_matches<<"\n";

      double stopTime = omp_get_wtime();
      double secsElapsed = stopTime - startTime; // that's all !
      std::cout<<"Matching time (no MM):"<<secsElapsed<<"\n";


      if (n_feat_frame_matches > min_matches_init_pose)
      {
        // Try to estimate the H and F from mathces

        // if ok...try more matches with model
        tracking_feat_cur_ref_ids = match_cur_ref_idx;


        std::cout<<"Init motion recovery...try next frame\n";
        // If fails try with next frame
        //selectNewFeaturesForTracking(ima, candidate_pts, candidate_pts_used, candidate_descs, mCurrentFrame->regions.get(), max_tracked_points - n_feat_frame_matches);

      }
      else
      {
        std::cout<<"Not enough matches for initialization process - setting this as REF frame\n";
        setReferenceSystemInitialization(mCurrentFrame);
      }

    }
  }

  size_t matchFramesFeatureMatchingNoMM
  (
    features::Regions * ref_feat_regions,
    features::PointFeatures & ref_feat_undist,
    features::Regions * candidate_feat_regions,
    features::PointFeatures & candidate_feat_undist,
    Hash_Map<size_t,size_t> & match_cur_ref_feat_ids,
    const size_t win_size = 50,
    const float ratio = 0.8
  )
  {
    std::cout<< "Match per frame!\n";

    // -------------------
    // -- Match reference-current candidates
    // -------------------
    std::cout<<"Matching candidates \n";

    std::vector<int>matches_ref_cur_idxs(ref_feat_undist.size(),-1);

    double startTime = omp_get_wtime();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t p_i=0; p_i < ref_feat_undist.size(); ++p_i)
    {
      //TODO: Get possible candidates through grid

      size_t best_idx = std::numeric_limits<size_t>::infinity();
      size_t second_best_idx = std::numeric_limits<size_t>::infinity();
      double best_distance = 30;//std::numeric_limits<double>::infinity();
      double second_best_distance = 30;//std::numeric_limits<double>::infinity();

      // Raw descriptors
      void * ref_pt_desc_raw;
      void * candidate_pt_desc_raw;

      // Get descriptor of reference point
      featureExtractor_->getDescriptorRaw(ref_feat_regions, p_i, &ref_pt_desc_raw);

      for (size_t c_i=0; c_i < candidate_feat_undist.size(); ++c_i)
      {
        // Check if points are close enough
        if ((candidate_feat_undist[c_i].coords() - ref_feat_undist[p_i].coords()).norm() > win_size)
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(candidate_feat_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance
        double distance = featureExtractor_->SquaredDescriptorDistance(ref_pt_desc_raw,candidate_pt_desc_raw);

        // Save if in best two
        if (distance < best_distance)
        {
          second_best_distance = best_distance;
          second_best_idx = best_idx;
          best_distance = distance;
          best_idx = c_i;
        }
        else if (distance < second_best_distance)
        {
          second_best_distance = distance;
          second_best_idx = c_i;
        }
      }

      // Detect best match
      if (best_idx != std::numeric_limits<size_t>::infinity())
      {
        if (second_best_idx != std::numeric_limits<size_t>::infinity())
        {
          if ((best_distance / second_best_distance) < ratio)
          {
            // Best is unique enough
            matches_ref_cur_idxs[p_i] = best_idx;
            //std::cout<<"II: P: "<<p_i<<": ("<<ref_feat_undist[p_i]<<") C: "<<best_idx<<": ("<<candidate_feat_undist[best_idx]<<") :: "<<best_distance<<" :: "<<(best_distance / second_best_distance)<<"\n";
          }
        }
        else
        {
          // Best is unique
          matches_ref_cur_idxs[p_i] = best_idx;
          //std::cout<<"III: P: "<<p_i<<": ("<<ref_feat_undist[p_i]<<") C: "<<best_idx<<": ("<<candidate_feat_undist[best_idx]<<") :: "<<best_distance<<" :: "<<(best_distance / second_best_distance)<<"\n";
        }
      }
    }

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Match candidates times:"<<secsElapsed<<"\n";


    // -------------------
    // -- Prune matches and insert
    // -------------------
    std::cout<<"Purging candidates and copy\n";

    startTime = omp_get_wtime();

    // Check that two are not matching the same point
    bool bOkMatch = true;
    for (size_t i=0; i<matches_ref_cur_idxs.size(); ++i)
    {
      if (matches_ref_cur_idxs[i] != -1)
      {
        for (size_t j=i+1; j<matches_ref_cur_idxs.size(); ++j)
        {
          // if value is doubled we delete both matches (mathces have to be unique)
          if (matches_ref_cur_idxs[i] == matches_ref_cur_idxs[j])
          {
            bOkMatch = false;
            matches_ref_cur_idxs[j] = -1;
          }
        }
        if (bOkMatch)
        {
          // Match
          match_cur_ref_feat_ids[matches_ref_cur_idxs[i]] = i;
        }
        else
        {
          // reset duplicate flag
          bOkMatch = true;

        }
      }
    }

    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Prune time:"<<secsElapsed<<"\n";

    return match_cur_ref_feat_ids.size();
  }

  bool detect
  (
    const image::Image<unsigned char> & ima,
    std::shared_ptr<Frame> &frame,
    const size_t min_count,
    const size_t max_count
  )
  {
    // Detect feature points
    size_t n_feats_detected = featureExtractor_->detect(ima,frame->regions,min_count,max_count);
    // Describe detected features
    featureExtractor_->describe(ima,frame->regions.get());
    // Undistort points
    frame->pts_undist = frame->regions->GetRegionsPositions();
    std::cout<<"PP: "<<frame->pts_undist.size();
    if (n_feats_detected > 0)
      return true;
    return false;
  }


};

} // namespace VO
} // namespace openMVG
