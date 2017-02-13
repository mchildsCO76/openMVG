
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

#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/tracking/Abstract_FeatureExtractor.hpp>
#include "nonFree/sift/SIFT_describer.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"

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

    max_dist_desc_d2 = 200*200;
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
    frame->regions_.reset(new RegionT);

    image_describer->Describe(ima, frame->regions_, nullptr);

//    RegionT * regionsCasted = dynamic_cast<RegionT*>(frame->regions_.get());
    /*
    RegionT * regionsCasted = dynamic_cast<RegionT*>(frame->regions_.get());

    const std::vector<unsigned char> scores = {30, 20, 10, 5};

    features::PointFeatures feats;
    // use a sequence of scores 'in order to deal with lighting change'
    for (const unsigned char fast_score : scores)
    {
      features::FastCornerDetector fastCornerDetector(9, fast_score);
      fastCornerDetector.detect(ima, feats);
      if (feats.size() > min_count)
      {
        // shuffle to avoid to sample only in one bucket
        if (max_count != 0)
        {
          std::random_shuffle(feats.begin(), feats.end());
          feats.resize(max_count); // cut the array to keep only a given count of features
        }
        regionsCasted->Features().resize(feats.size());
        std::copy(feats.begin(),feats.end(),regionsCasted->Features().begin());
        break;
      }
      else if (fast_score == scores[scores.size()-1] || feats.size() == min_count)
      {
        // If there is enough matches (or last detection step) we copy what we have
        regionsCasted->Features().resize(feats.size());
        std::copy(feats.begin(),feats.end(),regionsCasted->Features().begin());
        break;
      }
      feats.clear();
    }
*/
    // Estimate the uncertainty of each feature detection
    frame->pts_cov_.resize(frame->regions_->RegionCount());
    std::fill(frame->pts_cov_.begin(),frame->pts_cov_.end(), Eigen::Matrix2d::Identity());
    // Return number of detected features
    return frame->regions_->RegionCount();
  }

  bool describe
  (
    const image::Image<unsigned char> & ima,
    Frame * frame
  ) override
  {
    /*RegionT * regionsCasted = dynamic_cast<RegionT *>(frame->regions_.get());
    RegionT::FeatsT & pts = regionsCasted->Features();
    RegionT::DescsT & descs = regionsCasted->Descriptors();

    // Resize regions to fit new points
    descs.resize(regionsCasted->RegionCount());

    // Get descriptors of new points
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t i = 0; i < pts.size(); ++i)
    {
      const Vec2f & pt_pos = pts[i].coords();
      features::PickASDipole(ima, pt_pos(0), pt_pos(1), 10.5f, 0.0f, descs[i].data());
    }
    */
    return true;
  }

  void getDescriptorRaw
  (
    features::Regions * const regions,
    const size_t i,
    void ** desc
  ) override
  {
    //RegionT * regionsCasted = dynamic_cast<RegionT*>(regions);
    *desc = (dynamic_cast<RegionT *>(regions)->Descriptors()[i].data());
  }

  double SquaredDescriptorDistance
  (
      void * desc_A,
      void * desc_B
  ) override
  {

    openMVG::matching::L2_Vectorized<RegionT::DescriptorT::bin_type> metric;
    RegionT::DescriptorT::bin_type * d_A = static_cast<RegionT::DescriptorT::bin_type *>(desc_A);
    RegionT::DescriptorT::bin_type * d_B = static_cast<RegionT::DescriptorT::bin_type *>(desc_B);

    return metric(d_A, d_B, RegionT::DescriptorT::static_size);
  }

};

}
}
