
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

#include <openMVG/vsslam/detection/Abstract_FeatureExtractor.hpp>
#include <openMVG/vsslam/Frame.hpp>

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

    max_dist_desc_ = 200;
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
    const size_t min_count
  ) const override
  {
    // Cast region
    // Cast region
    std::unique_ptr<features::Regions> & frame_regions = frame->getRegionsRaw();
    frame_regions.reset(new RegionT);

    image_describer->Describe(ima, frame_regions, nullptr);

    const size_t n_features = frame_regions->RegionCount();
    // Estimate the uncertainty of each feature detection
    frame->pts_information_mat_.resize(n_features);
    // Scale of features
    frame->pts_scale_.resize(n_features);

    std::vector<RegionT::FeatureT> vec_features = dynamic_cast<RegionT *>(frame_regions.get())->Features();
    for (size_t i = 0; i < vec_features.size(); i++)
    {
      frame->pts_scale_[i] = vec_features[i].scale();
      frame->pts_information_mat_[i] = Eigen::Matrix2d::Identity() * (1.0 / 4);
    }

    // We assume sigma error of 1px
    //double sqrt_inv_sigma_px = 1; // (sqrt(1/1*1))
    //std::fill(frame->pts_information_mat_.begin(),frame->pts_information_mat_.end(), Eigen::Matrix2d::Identity()*sqrt_inv_sigma_px);
    // Return number of detected features
    return n_features;
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
