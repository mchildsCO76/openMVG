// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/types.hpp>
#include <software/VSSLAM/slam/Frame.hpp>

namespace openMVG  {
namespace VSSLAM  {

struct Abstract_FeatureExtractor
{
  virtual size_t detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count
  ) const = 0;

  virtual bool allocate
  (
    std::unique_ptr<features::Regions> & regions,
    size_t max_elements
  ) = 0;
  virtual bool resize
  (
    std::unique_ptr<features::Regions> & regions,
    size_t n_elements
  )=0;
  virtual bool describe
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & new_pts,
    std::unique_ptr<features::Regions> & regions
  )=0;
  virtual bool describe
  (
    const image::Image<unsigned char> & ima,
    features::PointFeature & pt,
    std::unique_ptr<void *> & desc
  )=0;
  virtual bool insert
  (
    std::unique_ptr<features::Regions> & regions,
    features::PointFeature & pt,
    std::unique_ptr<openMVG::Vecf> & desc,
    size_t idx
  )=0;
  virtual double SquaredDescriptorDistance
  (
    openMVG::Vecf * desc_A,
    openMVG::Vecf * desc_B
  )=0;
  virtual bool getDescriptorFromFrame
  (
    const std::shared_ptr<Frame> & ref_frame,
    const size_t feature_id,
    openMVG::Vecf * desc
  )=0;
/*
  virtual size_t getFrameMatching
  (
    std::shared_ptr<Frame> ref_frame,
    const image::Image<unsigned char> & current_ima,
    std::vector<features::PointFeature> & candidate_pts,
    const size_t window,
    const float ratio,
    std::unique_ptr<features::Regions> & new_feat_regions,
    Hash_Map<size_t,size_t> & feat_cur_new_matches_ids
  )=0;*/
};

} // namespace VSSLAM
} // namespace openMVG
