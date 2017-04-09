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
  float init_min_cos_angle_pt = 0.99995 ; // Min angle between rays for the point to be triangulated (0.99998 ~ 0.36deg; 0.99995 ~ 0.5deg;  0.9998 ~ 1.15deg)
  // Matching is 2D-2D
  float init_match_desc_ratio = 0.8; // Matching ration of descriptors
  float init_match_max_scale_ratio = 2; // Matching can be done only between points with max twice the scale difference


  size_t init_track_min_matches = 30; // Min number of points needed for initialization of tracking
  float default_reproj_thresh_sq = 4.0; // squared reprojection error used as default value for frames
  // ---------------
  // Matching settings
  // ---------------
  // Track with Motion Model
  float track_match_mm_desc_ratio = 0.9; // Matching ratio for matching using motion model
  float track_match_mm_max_scale_ratio = 1.2; // Matching can be done only between points with max twice the scale difference
  float track_mm_win_size = 5;

  // Track with Local Map
  float track_match_lm_desc_ratio = 0.8; // Matching ratio for matching using localMap
  float track_match_lm_max_scale_ratio = 2; // Matching can be done only between points with max twice the scale difference
  float track_lm_win_size = 5;
  float track_lm_pose_BA_inlier_ratio = 0.75; // If inlier ratio of matched map landmarks is below threshold we perform pose BA
  size_t track_local_map_size = 5;  // Number of best connected frames are used for constructing a local map and finding matches from map

  // Track with Reference frame
  float track_match_rf_desc_ratio = 0.7; // Matching ratio for matching using reference frame
  float track_rf_win_size = 10;
  float track_match_rf_max_scale_ratio = 2; // Matching can be done only between points with max twice the scale difference

  // Track with reference map
  float track_match_rm_desc_ratio = 0.7; // Matching ratio for matching using reference frame

  float track_match_rm_max_scale_ratio = std::numeric_limits<float>::infinity(); // Matching can be done only between points with max twice the scale difference


  // ---------------
  // Tracking settings
  // ---------------
  size_t track_mm_n_frames = 2; // Track points with mm from last n prev frames
  size_t track_min_matches = 10;  // Min matches needed for the tracking to be successful

  // ---------------
  // Triangulation settings
  // ---------------
  size_t triangulate_local_map_size = 10; // Number of best local frames used for triangulating the map
  size_t triangulate_min_obs_for_new_pt = 4; // Minimum observations required for new triangulation (only to start the process)
  float triangulate_match_lm_max_scale_ratio = 1.2; // Matching can be done only between points with max twice the scale difference
  float triangulate_match_epipolar_desc_ratio = 0.6;
  float track_epipolar_dist = 4;


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
  size_t last_relocalized_frame_id = 0;

  double chi2_px2_a005_thresh =  5.99; // chi2 (alpha=0.05,k=2)

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
    display_iterations = std::vector<size_t>(6,0);
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
    std::shared_ptr<Frame> current_frame,
    const image::Image<unsigned char> * mask = nullptr
  ) override;

  bool needNewKeyframe(Frame * frame,
      std::vector<std::unique_ptr<MapLandmark> > & vec_new_pts_3D);

  /// INITIALIZATION
  void resetTrackingInitialization();
  void clearTrackingInitializationData() {
    std::cout
        << "Tracker: [Initialization] Clear tracking initialization data!\n";
    if (init_ref_frame) {
      init_ref_frame.reset();
    }
  }

  void startTrackingInitialization(std::shared_ptr<Frame> & frame);

  // Returns true if current frame is used in the initialization process (also if its put as ref image)
  bool tryTrackingInitialization(const image::Image<unsigned char> & ima);

  /// TRACKING
  bool trackWithMotionModel();

  bool trackWithReferenceFrame();

  bool trackWithReferenceFrameLocalMap(Frame * frame_ref);

  bool trackLocalMap();

  void findNewLandmarks(Frame * frame, std::vector<std::unique_ptr<MapLandmark> > & vec_new_pts_3D);

  void markInliersInFrame(Frame * frame, Hash_Map<MapLandmark *,IndexT> & map_putative_matches_3D_pts_frame_idx,bool b_check_error, size_t association_id);
  double computeInlierRatioMatches(Frame * frame, Hash_Map<MapLandmark *,IndexT> & map_putative_matches_3D_pts_frame_idx);


  void removeOutliersInFrame(Frame * frame);

  //void removeOutliersInLocalFrames(std::vector<Frame *> & local_map_frames);

  void removeOutliersInNewTriangulatedPoints(Frame * frame, std::vector<std::unique_ptr<MapLandmark> > & vec_putative_new_pts_3D);


  // Detect and describe points in current frame
  // Goal is to get at least min_count features
  void detect
  (
    const image::Image<unsigned char> & ima,
    Frame * frame,
    const size_t & min_count,
    const image::Image<unsigned char> * mask = nullptr
  );
};

} // namespace VO
} // namespace openMVG
