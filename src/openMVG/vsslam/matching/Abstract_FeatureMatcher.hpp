// Copyright (c) 2017 Klemen ISTENIC

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/types.hpp>
#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>

namespace openMVG  {
namespace VSSLAM  {

class Abstract_FeatureMatcher
{
private:

public:
  virtual ~Abstract_FeatureMatcher(){};

  virtual void matching_AllAll_2D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    matching::IndMatches & vec_putative_matches_1_2_idx,
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) =0;
  virtual void matching_AllAll_3D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const std::vector<MapLandmark *> & vec_3D_pts,
    const Frame * frame,
    Hash_Map<MapLandmark* ,IndexT> & matches_3D_pts_frame_idx,
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) =0;
  virtual void matching_Projection_3D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame,
    const std::vector<MapLandmark *> & vec_3D_pts,
    Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_idx,
    const size_t win_size = 15,
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) =0;
  virtual void matching_Projection_3D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_idx,
    const size_t max_px_d = 15, //20*20
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) =0;
  virtual void matching_EpipolarLine_2D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    const Mat3 & F_21,
    Hash_Map<IndexT,IndexT> & map_matches_1_2_idx,
    const size_t max_px_epi_d = 4,
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) =0;

  // Return radius depending on the angle between average viewing angle of a point and current viewing angle of observation
  float radiusByViewingAngle(const float cos_view_angle)
  {
    // arccos(0.998) ~ 3.6deg arccos(0.9998) ~ 1.15 deg
    if (cos_view_angle > 0.9995)
    {
      return 1.0;
    }
    else
    {
      //return 1.0 + (1-cos_view_angle)*4.0;
      return 2.0;
    }
  }

};

}
}
