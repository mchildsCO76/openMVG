
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/types.hpp>
#include <numeric>

#include <memory>
#include "openMVG/features/features.hpp"
#include "openMVG/features/image_describer.hpp"
#include "openMVG/image/image.hpp"

#include "nonFree/sift/SIFT_describer.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"

#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/detection/Abstract_FeatureExtractor.hpp>

namespace openMVG  {
namespace VSSLAM  {

class Feat_Extractor_SIFT : public Abstract_FeatureExtractor
{
  using RegionT = features::SIFT_Regions;

  std::unique_ptr<features::Image_describer> image_describer;
public:
  Feat_Extractor_SIFT(): Abstract_FeatureExtractor()
  {
    //image_describer.reset(new features::SIFT_Image_describer
    //  (features::SIFT_Image_describer::Params(), false));

    image_describer.reset(
      new features::SIFT_Anatomy_Image_describer(features::SIFT_Anatomy_Image_describer::Params()));

    max_dist_desc_ = 100;
    //max_dist_desc_ = std::numeric_limits<float>::infinity();
  }

  size_t getDescriptorLength() const override
  {
    return RegionT::DescriptorT::static_size;
  }


  // suggest new feature point for tracking (count point are kept)
  size_t detect
  (
    const image::Image<unsigned char> & ima,
    Frame * frame,
    const size_t min_count,
    const size_t max_count
  ) const override
  {
    // Cast region
    // Cast region
    std::unique_ptr<features::Regions> & frame_regions = frame->getRegionsRaw();
    frame_regions.reset(new RegionT);

    image_describer->Describe(ima, frame_regions, nullptr);

    // Estimate the uncertainty of each feature detection
    frame->pts_cov_.resize(frame_regions->RegionCount());
    std::fill(frame->pts_cov_.begin(),frame->pts_cov_.end(), Eigen::Matrix2d::Identity());
    // Return number of detected features
    return frame_regions->RegionCount();
  }

  bool describe
  (
    const image::Image<unsigned char> & ima,
    const Frame * frame
  ) const override
  {
    return true;
  }


  void getDescriptorRaw
  (
    features::Regions * const regions,
    const IndexT i,
    void ** desc
  ) const override
  {
    //RegionT * regionsCasted = dynamic_cast<RegionT*>(regions);
    *desc = (dynamic_cast<RegionT *>(regions)->Descriptors()[i].data());
  }

  double SquaredDescriptorDistance
  (
      void * desc_A,
      void * desc_B
  ) const override
  {

    openMVG::matching::L2_Vectorized<RegionT::DescriptorT::bin_type> metric;
    RegionT::DescriptorT::bin_type * d_A = static_cast<RegionT::DescriptorT::bin_type *>(desc_A);
    RegionT::DescriptorT::bin_type * d_B = static_cast<RegionT::DescriptorT::bin_type *>(desc_B);

    return metric(d_A, d_B, RegionT::DescriptorT::static_size);
  }

};

}
}
