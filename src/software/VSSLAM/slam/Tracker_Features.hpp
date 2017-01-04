// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>

#include <software/VSSLAM/slam/Frame.hpp>
#include <software/VSSLAM/slam/Abstract_Tracker.hpp>
#include <software/VSSLAM/slam/Abstract_FeatureExtractor.hpp>

namespace openMVG  {
namespace VSSLAM  {

struct Tracker_Features : public Abstract_Tracker
{
  /// Feature extractor
  std::unique_ptr<Abstract_FeatureExtractor> featureExtractor_;

  // Initialization
  std::shared_ptr<Frame> init_ref_frame;


  void setFeatureExtractor(std::shared_ptr<Abstract_FeatureExtractor> & featExtractor)
  {
    featureExtractor_.reset(featExtractor.get());
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
    // Allocate space for all features
    featureExtractor_->allocate(mCurrentFrame,max_tracked_points);

    std::cout<<"Tracks from prev frame: "<<prev_tracked_ids.size()<<"\n";

    // Check if we have tracks from before
    if (!prev_tracked_ids.empty())
    {
      std::cout<<"TRY TRACKING!\n";

      // Try to detect features from previous point
      std::vector<features::PointFeature> new_pts;
      // Detect features - all possible
      bool bDetect = detect(ima,new_pts,0);
      std::cout<<"Detected features: "<<new_pts.size()<<"\n";
      // if successful all are tracks
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
          std::cout<<"No motion model tracking!\n";
          // We do not have any motion model - normal search around points
          size_t win_size = 100;



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
          if (trackingStatus == TRACKING_STATUS::INIT)
          {
            resetSystemInitialization();
            systemInitialization();
          }
              // systemInitialization();
              // set this as new init ref
              // trackingStatus = TRACKING_STATUS::INIT;
        }
      }
      else
      {
        prev_tracked_ids.clear();
        if (trackingStatus == TRACKING_STATUS::INIT)
        {
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
      // No tracks in prev frames -> we detect new tracks
      std::vector<features::PointFeature> new_pts;
      // We detect enough new points to fill to max_tracked_points
      bool bDetect = detect(ima,new_pts,max_tracked_points);
      std::cout<<"N of proposed features: "<<new_pts.size()<<"\n";

      // if successful all are tracks
      if (bDetect)
      {
        std::unique_ptr<features::Regions> regions;
        // Retrieve descriptors and add to frame
        featureExtractor_->describe(ima,new_pts,regions);
        // Insert new features into frame
        featureExtractor_->insert(mCurrentFrame,regions);
        // Set tracks as tracked
        prev_tracked_ids.resize(mCurrentFrame->regions->RegionCount());
        std::iota(prev_tracked_ids.begin(), prev_tracked_ids.end(), 0);
        // Set tracker to lost (has points to track but its not actually tracking yet)
        std::cout<<"Feat after insertion "<<mCurrentFrame->regions->RegionCount()<<" \n";

        // If system is not initialized we signal that we can start initialization or else just lost procedure
        if (trackingStatus == TRACKING_STATUS::NOT_INIT)
        {
          trackingStatus = TRACKING_STATUS::INIT;
          systemInitialization();
        }
        else
        {
          trackingStatus = TRACKING_STATUS::LOST;
        }
      }
      else
      {
        // Unsuccessful detection of features
        //prev_tracked_ids.clear();
        std::cout<<"Feat FAILED \n";
        // If system was previously initialized we just put it to idle
        // TODO: check if its idle too long it has to become not initialized
        if (trackingStatus == TRACKING_STATUS::NOT_INIT)
        {
          trackingStatus = TRACKING_STATUS::NOT_INIT;
        }
        else if (trackingStatus == TRACKING_STATUS::INIT)
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

    std::cout<<"Track STATUS: ";
    switch (trackingStatus)
    {
    case Abstract_Tracker::TRACKING_STATUS::OK:
      std::cout<<"OK\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::LOST:
      std::cout<<"LOST\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::INIT:
      std::cout<<"INIT\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::NOT_INIT:
      std::cout<<"NOT_INIT\n";
      break;
    case Abstract_Tracker::TRACKING_STATUS::IDLE:
      std::cout<<"IDLE\n";
      break;
    }

    std::cout<<"Tracks available: "<<prev_tracked_ids.size()<<"\n";

    // Update frame structure in tracker
    mPrevFrame.swap(mCurrentFrame);
    mCurrentFrame.reset();

    // Return if tracking is ok
    return trackingStatus==TRACKING_STATUS::OK;
  }

  // suggest new feature point for tracking (count point are kept)
  bool detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count
  ) const override
  {
    return featureExtractor_->detect(ima,pt_to_track,count);
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
  // Try to initialize system
  void systemInitialization
  (
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

        // if successful
          // trackingStatus = TRACKING_STATUS::OK;
        // if not

        //return true;
      }
    }
  }



};

} // namespace VO
} // namespace openMVG
