// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/types.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"

#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"
#include "openMVG/matching/regions_matcher.hpp"


#include <openMVG/vsslam/system/Frame.hpp>
#include <openMVG/vsslam/vsslam_parameters.hpp>
#include <openMVG/vsslam/features/Abstract_Feature_Matcher.hpp>

namespace openMVG {
namespace vsslam {

class Feat_Matcher_Regions : public Abstract_Feature_Matcher
{
private:
  matching::EMatcherType e_matcher_type_;

public:
   Feat_Matcher_Regions
   (
     std::shared_ptr<VSSLAM_Parameters> & params,
     Abstract_Feature_Extractor * feat_extractor,
     matching::EMatcherType matcher_type = matching::EMatcherType::BRUTE_FORCE_L2
   ) : Abstract_Feature_Matcher(params)
   {
     e_matcher_type_ = matcher_type;
   }

  void matching_AllAll_2D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor,
    const Frame * frame_1,
    const Frame * frame_2,
    matching::IndMatches & vec_putative_matches_1_2_idx,
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    matching::Matcher_Regions_Database matcher(e_matcher_type_, *frame_1->getRegions());
    if (frame_1->getRegions()->Type_id() != frame_2->getRegions()->Type_id())
    {
      return;
    }

    //matching::IndMatches vec_putatives_matches;
    matcher.Match(f_desc_ratio, *frame_2->getRegions(), vec_putative_matches_1_2_idx);
  }

  void matching_AllAll_3D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor,
    const std::vector<MapLandmark *> & vec_landmarks,
    const Frame * frame,
    Hash_Map<MapLandmark* ,IndexT> & matches_landmark_frame_idx,
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {

  }
  void matching_Projection_3D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    Hash_Map<MapLandmark *,IndexT> & matches_frame_1_2_idx,
    const float f_max_px_d = 15, //20*20
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {

  }

  void matching_Projection_3D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor,
    const std::vector<MapLandmark *> & vec_landmarks,
    const Frame * frame,
    Hash_Map<MapLandmark* ,IndexT> & matches_landmark_frame_idx,
    const float f_max_px_d = 15,
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {

  }


  void matching_Epipolar_2D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    const Mat3 & F_21,
    Hash_Map<IndexT,IndexT> & map_matches_1_2_idx,
    const float f_max_px_d = 4,
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {

  }
};

}
}
