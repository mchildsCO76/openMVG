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

  //
  size_t min_init_ref_tracks = 30; // Min number of tracks reference frame for initialization has to have
  size_t min_frame_tracks = 10; // Min number of tracks frame has to have
  size_t max_frame_tracks = 0; // Max number of feats detected in a frame (0 - unlimited)
  size_t min_matches_init_pose = 30; // Min number of matches for init pose estimation

  // Tracking possibilities
  std::vector<features::PointFeature> candidate_pts;
  std::vector<bool> candidate_pts_used;


  Tracker_Features(){}

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
    bool bNewReferenceframe = true;

    feat_cur_prev_matches_ids.clear();
    // Set current frame
    mCurrentFrame = current_frame->share_ptr();
    // Allocate space for all features
    featureExtractor_->allocate(mCurrentFrame->regions,max_tracked_points);


    if (trackingStatus == TRACKING_STATUS::NOT_INIT)
    {
      // Detect new features and add frame as init reference frame
      std::cout<<"DETECT FROM STRACH A!\n";
      // We detect enough new points to fill to max_tracked_points
      bool bDetect = detect(ima,candidate_pts,max_tracked_points,min_init_ref_tracks);
      candidate_pts_used.resize(candidate_pts.size());
      std::fill(candidate_pts_used.begin(),candidate_pts_used.end(),false);

      std::cout<<"Candidate features: "<<candidate_pts.size()<<" :: "<<candidate_pts_used.size()<<"\n";

      // if successful we add them as features in the current frame
      if (bDetect)
      {
        // Retrieve descriptors and add to frame
        featureExtractor_->describe(ima,candidate_pts,mCurrentFrame->regions);
        std::cout<<"Feats after insertion: "<<mCurrentFrame->getTracksSize()<<" \n";

        trackingStatus = TRACKING_STATUS::INIT;
        setReferenceSystemInitialization(mCurrentFrame);

      }
      else
      {
        // Unsuccessful detection of features
        std::cout<<"Feat detection FAILED \n";
        // We go to next frame
      }
    }
    else if (trackingStatus == TRACKING_STATUS::INIT)
    {
      std::cout<<"TRY TO TRACK A!\n";

      // Detect features from frame and try to match to reference initialization frame
      // Detect all possible (hence 0)
      bool bDetect = detect(ima,candidate_pts,max_frame_tracks,min_frame_tracks);
      candidate_pts_used.resize(candidate_pts.size());
      std::fill(candidate_pts_used.begin(),candidate_pts_used.end(),false);

      std::cout<<"Candidate features: "<<candidate_pts.size()<<" :: "<<candidate_pts_used.size()<<"\n";
      std::cout<<"Init ref features: "<<init_ref_frame->getTracksSize()<<"\n";

      // if we detect enough we try to match to ref init frame
      if (bDetect)
      {
        trySystemInitialization(ima);

      }
      else
      {
        // not enough features detected
        resetSystemInitialization();
      }


    }

/*


    // Get number of tracks from prev frame
    size_t n_tracks_lastFrame = 0;
    if (mPrevFrame)
    {
      n_tracks_lastFrame = mPrevFrame->getTracksSize();
    }
    std::cout<<"Tracks from prev frame: "<<n_tracks_lastFrame<<"\n";



    // Check if we have tracks from before
    if (n_tracks_lastFrame!=0)
    {
      std::cout<<"TRY TRACKING!\n";

      // Detect all possible (hence 0)
      bool bDetect = detect(ima,candidate_pts,0);

      std::cout<<"Detected features: "<<candidate_pts.size()<<"\n";

      // if we detect enough we try to match to last frame
      if (bDetect)
      {
        // Check if we can use motion model
        if (trackingStatus == TRACKING_STATUS::OK)
        {
          // We use the model to predict the position of the tracks
          std::cout<<"Motion model tracking!\n";

          // if tracking sucessful
          trackingStatus = TRACKING_STATUS::OK;
          // else trackingStatus = TRACKING_STATUS::LOST;

        }
        else
        {

          if (trackingStatus == TRACKING_STATUS::INIT)
          {
            systemInitialization(&ima,&candidate_pts);
          }
          else
          {

          }




          // We do not have any motion model - normal search for matches from reference frame
          // (in window around last known position)
          // Search parameters
         /* size_t win_size = 30;
          float desc_ratio = 0.8;
          std::cout<<"No motion model tracking!\n";

          // Map between current(matched_regions) and ref_regions ids
          Hash_Map<size_t,size_t> &feat_cur_ref_matches_ids;

          // Points of new frame
          std::unique_ptr<features::Regions> matched_regions;

          // Find matches between current and last reference frame
          size_t n_feat_frame_matches = getMatchesFrameFeatureMatchingNoMotionModel(mPrevFrame,candidate_pts,win_size,desc_ratio,new_regions,feat_frame_matches_ids);
          std::cout<<"Candidate frame matches: "<<n_feat_frame_matches<<"\n";

          if (n_frame_matches < 100)
          {


          }
          else
          {
            if (trackingStatus == TRACKING_STATUS::INIT)
            {
              resetSystemInitialization();
              systemInitialization();
            }
            else
            {
              trackingStatus = TRACKING_STATUS::LOST;
            }
          }*/
          // if successful
            // if lost
              // trackingStatus = TRACKING_STATUS::OK;
            // if init
              // systemInitialization();
              // try to initialize
              // trackingStatus = TRACKING_STATUS::INIT;
          // if not successful
            // if lost
              // trackingStatus = TRACKING_STATUS::LOST;
            // if init
              // systemInitialization();
              // set this as new init ref
              // trackingStatus = TRACKING_STATUS::INIT;
   /*     }
      }
      else
      {
        std::cout<<"Feat detection FAILED \n";
        // We cant detect features
        if (trackingStatus == TRACKING_STATUS::INIT)
        {
          resetSystemInitialization();
          trackingStatus = TRACKING_STATUS::NOT_INIT;
        }
        else
        {
          trackingStatus = TRACKING_STATUS::IDLE;
        }
      }
    }
    else
    {
      std::cout<<"DETECT FROM STRACH!\n";
      // We detect enough new points to fill to max_tracked_points
      bool bDetect = detect(ima,candidate_pts,max_tracked_points);
      std::cout<<"Candidate features: "<<candidate_pts.size()<<"\n";

      // if successful we add them as features in the current frame
      if (bDetect)
      {
        // Retrieve descriptors and add to frame
        featureExtractor_->describe(ima,candidate_pts,mCurrentFrame->regions);
        // Insert new features into frame
        //featureExtractor_->insert(mCurrentFrame,regions);
        // Set tracks as tracked
        //prev_tracked_ids.resize(mCurrentFrame->regions->RegionCount());
        //std::iota(prev_tracked_ids.begin(), prev_tracked_ids.end(), 0);
        // Set tracker to lost (has points to track but its not actually tracking yet)
        std::cout<<"Feats after insertion: "<<mCurrentFrame->regions->RegionCount()<<" \n";

        // If system is not initialized we signal that we can start initialization or else just lost procedure
        if (trackingStatus == TRACKING_STATUS::NOT_INIT)
        {
          trackingStatus = TRACKING_STATUS::INIT;
          systemInitialization();
        }
        else
        {
          trackingStatus = TRACKING_STATUS::LOST; // we have features but no previous ones -> need relocalization
        }
      }
      else
      {
        // Unsuccessful detection of features
        //prev_tracked_ids.clear();
        std::cout<<"Feat detection FAILED \n";
        // TODO: check if its idle too long it has to become not initialized
        if (trackingStatus == TRACKING_STATUS::NOT_INIT)
        {
          // we keep non_init status
          //trackingStatus = TRACKING_STATUS::NOT_INIT;
        }
        else if (trackingStatus == TRACKING_STATUS::INIT)
        {
          // If we are in the initialization process and get empty frame
          // Reset and start again from beginning with next frame
          resetSystemInitialization();
          trackingStatus = TRACKING_STATUS::NOT_INIT;
        }
        else
        {
          // If system was previously initialized we just put it to idle - no features in this frame but it can just be temporaty
          trackingStatus = TRACKING_STATUS::IDLE;
        }
      }
    }
      */

    printTrackingStatus();

    std::cout<<"Tracks available: "<<mCurrentFrame->regions->RegionCount()<<"\n";

    // Update frame structure in tracker
    //if (bNewReferenceframe)
    //{
      mPrevFrame.swap(mCurrentFrame);
    //}

    mCurrentFrame.reset();
    // Candidate features in current frame
    candidate_pts.clear();
    candidate_pts_used.clear();

    // Return if tracking is ok
    return trackingStatus==TRACKING_STATUS::OK;
  }

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
    // Check if we have enough features in the frame
    if (trackingStatus == Abstract_Tracker::TRACKING_STATUS::INIT)
    {
      // Set current frame as the new reference frame for initialization
      std::cout<<"Set REF initialization frame\n";
      init_ref_frame = frame->share_ptr();
    }
  }

  void trySystemInitialization
  (
    const image::Image<unsigned char> & ima
  )
  {
    std::cout<<"System initialized process B!\n";
    // Check if we have enough features in the frame
    if (trackingStatus == Abstract_Tracker::TRACKING_STATUS::INIT)
    {
      std::cout<<"Try to match with ref init frame\n";

      // Matching with the init reference frame
      size_t win_size = 50;
      float desc_ratio = 0.8;

      // Map between current(matched_regions) and ref_regions ids
      //Hash_Map<size_t,size_t> feat_cur_ref_matches_ids;
      size_t n_feat_frame_matches = matchFramesFeatureMatchingNoMM(init_ref_frame,ima,candidate_pts,candidate_pts_used,win_size,desc_ratio,mCurrentFrame->regions,feat_cur_prev_matches_ids);

      if (n_feat_frame_matches > min_matches_init_pose)
      {
        // Try to estimate the H and F

        // If fails try with next frame
      }
      else
      {
        // Not enough matches
        // Select new features for the frame and add as initialization reference frame


        std::cout<<"Not enough matches for initialization process - setting this as REF frame\n";
        setReferenceSystemInitialization(mCurrentFrame);
      }
    }
  }

  /*
  // Try to initialize system
  void systemInitialization
  (
    const image::Image<unsigned char> * ima = nullptr,
    std::vector<features::PointFeature> * candidate_pts = nullptr
  )
  {
    std::cout<<"System initialized process!\n";

    // Check if we have enough features in the frame
    if (trackingStatus == Abstract_Tracker::TRACKING_STATUS::INIT)
    {
      // Check if we have a reference init frame
      if (!init_ref_frame)
      {
        // Set current frame as the new reference frame for initialization
        std::cout<<"Set REF initialization frame\n";
        init_ref_frame = mCurrentFrame->share_ptr();
      }
      else
      {
        // Tracking from previous frame was successful
        std::cout<<"Try to initialize\n";

        // Matching with the init reference frame
        size_t win_size = 10;
        float desc_ratio = 0.8;

        // Map between current(matched_regions) and ref_regions ids
        Hash_Map<size_t,size_t> feat_cur_ref_matches_ids;

        // Search with no motion model
        size_t n_feat_frame_matches = getMatchesFrameFeatureMatchingNoMotionModel(init_ref_frame,*ima,*candidate_pts,win_size,desc_ratio,mCurrentFrame->regions,feat_cur_ref_matches_ids);

        // Select new tracks amongs the rest
        std::unique_ptr<features::Regions> new_feat_regions;
        std::vector<features::PointFeature> new_track_pts;
        selectNewFeaturesForTracking(*candidate_pts,feat_cur_ref_matches_ids,new_track_pts,(max_tracked_points - n_feat_frame_matches));
        std::cout<<"New features detected: "<<new_track_pts.size()<<"\n";

        if (new_track_pts.size() > 0)
        {
          // Get descriptors of new points
          featureExtractor_->describe(*ima,new_track_pts,new_feat_regions);
          std::cout<<"New features described: "<<new_feat_regions->RegionCount()<<"\n";
          featureExtractor_->insert(mCurrentFrame,new_feat_regions);
          std::cout<<"New features of frame: "<<new_feat_regions->RegionCount()<<"\n";
        }

        std::cout<<"Candidate frame matches: "<<n_feat_frame_matches<<"\n";
        std::cout<<"Candidate frame features: "<<mCurrentFrame->regions->RegionCount()<<"\n";

        if (n_feat_frame_matches > 100)
        {
          // try estimate pose
        }
        else
        {
          // not enough matches
          // make this a new reference frame
          resetSystemInitialization();
          systemInitialization();
        }
        // if successful
          // trackingStatus = TRACKING_STATUS::OK;
        // if not

        //return true;
      }
    }
  }
*/

  size_t matchFramesFeatureMatchingNoMM
  (
    std::shared_ptr<Frame> & ref_frame,
    const image::Image<unsigned char> & current_ima,
    std::vector<features::PointFeature> & candidate_pts,
    std::vector<bool> & candidate_pts_used,
    size_t win_size,
    float ratio,
    std::unique_ptr<features::Regions> & new_feat_regions,
    Hash_Map<size_t,size_t> & feat_cur_ref_matches_ids
  )
  {
    // For each point in frame find candidates that are close enough by distance

      // Set of id of points from ref and all possible match_id of features from current frame
      std::vector<std::set<size_t> > vec_ref_cur_candidates(ref_frame->getTracksSize());

      // Get pos of features in the reference frame
      features::PointFeatures ref_points_pos = ref_frame->regions->GetRegionsPositions();

      std::cout<<"Finding possible candidate matches (ref->cur) "<<ref_frame->frameId_<<"\n";
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for (size_t p_i=0; p_i < ref_points_pos.size(); ++p_i)
      {
        for (size_t c_i=0; c_i < candidate_pts.size(); ++c_i)
        {
          if ((candidate_pts[c_i].coords() - ref_points_pos[p_i].coords()).norm() < win_size)
          {
            vec_ref_cur_candidates[p_i].emplace(c_i);
            candidate_pts_used[c_i] = true;
          }
        }
      }
      // Compute how many new points we will investigate
      int n_to_check_newpts = std::count(candidate_pts_used.begin(), candidate_pts_used.end(), true);
      std::cout<<"Need to check: "<<n_to_check_newpts<<" new points\n";



      // Get only useful candidates
      std::cout<<"Describe useful candidates\n";
      std::unique_ptr<features::Regions> candidate_regions;

      featureExtractor_->describe(current_ima,candidate_pts,candidate_regions);


      // Try to match point in reference frame with the identified candidates
      std::vector<int>matches_ref_cur_idxs(ref_points_pos.size(),-1);

      std::cout<<"Matching candidates "<<candidate_regions->RegionCount()<<"\n";
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for (size_t p_i=0; p_i < ref_points_pos.size(); ++p_i)
      {
        if (vec_ref_cur_candidates[p_i].empty())
          continue;

        size_t best_idx = std::numeric_limits<size_t>::infinity();
        size_t second_best_idx = std::numeric_limits<size_t>::infinity();
        double best_distance = 30;//std::numeric_limits<double>::infinity();
        double second_best_distance = 30;//std::numeric_limits<double>::infinity();

        for (size_t c_i : vec_ref_cur_candidates[p_i])
        {
          double distance = ref_frame->regions->SquaredDescriptorDistance(p_i,candidate_regions.get(),c_i);
          //double distance = featureExtractor_->SquaredDescriptorDistance(ref_pt_desc,candidate_descs[c_i].get());

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
        if (best_idx != std::numeric_limits<size_t>::infinity())
        {
          if (second_best_idx != std::numeric_limits<size_t>::infinity())
          {
            if ((best_distance / second_best_distance) > ratio)
            {
              // Best is unique enough
              matches_ref_cur_idxs[p_i] = best_idx;
            }
          }
          else
          {
            // Best is unique
            matches_ref_cur_idxs[p_i] = best_idx;
          }
        }
      }


      std::cout<<"Purging candidates\n";
      size_t n_matches=0;
      // Check that two are not matching the same point
      for (size_t i=0; i<ref_points_pos.size(); ++i)
      {
        if (matches_ref_cur_idxs[i] != -1)
        {
          bool bDuplicate = false;
          for (size_t j=i+1; j<ref_points_pos.size(); ++j)
          {
            // if value is doubled we delete both matches (mathces have to be unique)
            if (matches_ref_cur_idxs[i] == matches_ref_cur_idxs[j])
            {
              bDuplicate = true;
              matches_ref_cur_idxs[j] = -1;
            }
          }
          if (!bDuplicate)
          {
            //Unique match -> add to the features
            feat_cur_ref_matches_ids[n_matches] = i;
            n_matches++;
          }
          else
          {
            matches_ref_cur_idxs[i] = -1;
          }
        }
      }

      std::cout<<"Copy feature data to frame "<<n_matches<<"\n";
      //featureExtractor_->resize(new_feat_regions,n_matches);

      std::cout<<"Copy feature data to frame "<<n_matches<<"\n";
      //#ifdef OPENMVG_USE_OPENMP
      //#pragma omp parallel for schedule(dynamic)
      //#endif
      for (size_t i=0; i<n_matches; ++i)
      {
        const size_t & ref_pt_idx = feat_cur_ref_matches_ids[i];
        const size_t & candidate_pt_idx = matches_ref_cur_idxs[ref_pt_idx];
        //featureExtractor_->insert(new_feat_regions,candidate_pts[candidate_pt_idx],candidate_descs[candidate_pt_idx],i);

        candidate_regions->CopyRegion(candidate_pt_idx,new_feat_regions.get());
      }

      std::cout<<"Features copied: " <<new_feat_regions->RegionCount()<< " :: "<<n_matches<<"\n";






    // Get feature candidates
   // size_t n_frame_matches = featureExtractor_->getFrameMatching(frame,current_ima,candidate_pts,win_size,ratio,new_feat_regions,feat_cur_ref_matches_ids);

    return n_matches;
  }

/*
  size_t getMatchesFrameFeatureMatchingNoMotionModel
  (
    std::shared_ptr<Frame> & frame,
    const image::Image<unsigned char> & current_ima,
    std::vector<features::PointFeature> & candidate_pts,
    size_t win_size,
    float ratio,
    std::unique_ptr<features::Regions> & new_feat_regions,
    Hash_Map<size_t,size_t> & feat_cur_ref_matches_ids
  )
  {
    // Get feature candidates
    size_t n_frame_matches = featureExtractor_->getFrameMatching(frame,current_ima,candidate_pts,win_size,ratio,new_feat_regions,feat_cur_ref_matches_ids);

    return n_frame_matches;
  }*/
/*
  void selectNewFeaturesForTracking
  (
    const std::vector<features::PointFeature> & candidate_pts,
    const Hash_Map<size_t,size_t> & feat_frame_matches_ids,
    std::vector<features::PointFeature> & new_track_pts,
    size_t n_new_tracks
  )
  {
    // Get random order
    std::vector<int> possible_ids(new_pts.size());
    std::iota(possible_ids.begin(), possible_ids.end(), 0);
    std::random_shuffle(possible_ids.begin(), possible_ids.end());
    // Reserve space for new ids
    new_track_pts.reserve(n_new_tracks);

    for (size_t pos_id: possible_ids)
    {
      // Check if it has not been already matched
      if (feat_frame_matches_ids.count(pos_id) == 0)
      {
        features::PointFeature pf(new_pts[pos_id]);
        new_track_pts.push_back(pf);
        n_new_tracks--;
        if (n_new_tracks==0)
        {
          break;
        }
      }
    }
  }
*/

};

} // namespace VO
} // namespace openMVG
