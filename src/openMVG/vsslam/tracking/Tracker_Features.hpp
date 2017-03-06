// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>
#include <tuple>

#include <openMVG/types.hpp>
#include "openMVG/multiview/triangulation.hpp"

#include <openMVG/sfm/sfm.hpp>

#include <openMVG/vsslam/tracking/Abstract_Tracker.hpp>
#include <openMVG/vsslam/detection/Abstract_FeatureExtractor.hpp>
#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/matching/Abstract_FeatureMatcher.hpp>
#include <openMVG/vsslam/tracking/PoseEstimation.hpp>

namespace openMVG  {
namespace VSSLAM  {

class Tracker_Features : public Abstract_Tracker
{
private:
  /// Feature extractor
  Abstract_FeatureExtractor * featureExtractor_;
  /// Feature matcher
  Abstract_FeatureMatcher * featureMatcher_;

  // ---------------
  // Parameters
  // ---------------
  size_t max_tracked_points = 1500;

  // ---------------
  // Tracking init settings
  // ---------------
  size_t init_track_min_matches = 30; // Min number of points needed for initialization of tracking
  float init_min_cos_angle_pt = 0.99998; // Min angle between rays for the point to be triangulated (0.99998 ~ 0.36deg; 0.99995 ~ 0.5deg;  0.9998 ~ 1.15deg)

  float default_reproj_thresh_2 = 16.0; // squared reprojection error used as default value for frames
  // ---------------
  // Matching settings
  // ---------------
  float init_match_desc_ratio = 0.8; // Matching ration of descriptors
  float track_match_desc_ratio = 0.8;
  float track_epipolar_desc_ratio = 0.6;

  float track_mm_win_size = 15;
  float track_rf_win_size = 15;
  float track_local_map_win_size = 10;
  float track_epipolar_dist = 4;

  // ---------------
  // Tracking settings
  // ---------------
  size_t track_min_matches = 40;  // Min matches needed for the tracking to be successful
  size_t track_local_map_size = 5;  // Number of best connected frames are used for constructing a local map and finding matches from map
  size_t triangule_local_map_size = 2; // Number of best local frames used for triangulating the map

  // ---------------
  // Mapping settings
  // ---------------
  // Set  in cartographer.hpp
  /*
  size_t init_min_map_pts = 30; // Min points needed for successful map initialization
  size_t min_obs_per_landmark = 3;  // Min observations needed for the landmark to be valid
  size_t min_landmark_per_frame = 3;  // Min landmarks per frame to consider frame defined
  size_t max_frames_inactive_local_landmark = 5;  // Remove local landmarks after they are not observed for X frames
   */


  double vec_times[10];

public:
  // ---------------
  // Parameters
  // ---------------
  Tracker_Features
  (
    Abstract_FeatureExtractor * featExtractor,
    Abstract_FeatureMatcher * featMatcher,
    const size_t max_features_tracked = 1500
  ) :
    featureExtractor_(featExtractor),
    featureMatcher_(featMatcher),
    max_tracked_points(max_features_tracked)
  {
  }

  ~Tracker_Features() = default;

  void setFeatureExtractor(Abstract_FeatureExtractor * featExtractor)
  {
    featureExtractor_ = featExtractor;
  }

  void setFeatureMatcher(Abstract_FeatureMatcher * featMatcher)
  {
    featureMatcher_ = featMatcher;
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
  ) override;



  /// INITIALIZATION
  void resetSystemInitialization();

  void setReferenceSystemInitialization(std::shared_ptr<Frame> & frame);

  // Returns true if current frame is used in the initialization process (also if its put as ref image)
  bool trySystemInitialization(const image::Image<unsigned char> & ima);

  /// TRACKING
  bool trackWithMotionModel();

  bool trackWithReferenceFrame();

  void trackLocalMap();

  void triangulateNewLandmarks(std::vector<std::unique_ptr<MapLandmark> > & vec_new_pts_3D);

  void checkReprojectionAndMarkInliersInFrame(Frame * frame, Hash_Map<MapLandmark *,IndexT> & map_putative_matches_3D_pts_frame_idx);

  void markInliersInFrame(Frame * frame, Hash_Map<MapLandmark *,IndexT> & map_putative_matches_3D_pts_frame_idx);

  //void removeOutliersInLocalFrames(std::vector<Frame *> & local_map_frames);

  void removeOutliersInNewTriangulatedPoints(std::vector<std::unique_ptr<MapLandmark> > & vec_putative_new_pts_3D);


  // Detect and describe points in current frame
  bool detect
  (
    const image::Image<unsigned char> & ima,
    Frame * frame,
    const size_t & min_count,
    const size_t & max_count
  );
};

} // namespace VO
} // namespace openMVG
