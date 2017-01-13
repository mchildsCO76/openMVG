
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/types.hpp>
#include <numeric>

#include <software/VSSLAM/slam/Frame.hpp>
#include <software/VSSLAM/slam/Abstract_FeatureExtractor.hpp>


namespace openMVG  {
namespace VSSLAM  {

struct Feat_Extractor_FastDipole : public Abstract_FeatureExtractor
{
  using RegionT = features::FAST_Dipole_Regions;

  // suggest new feature point for tracking (count point are kept)
  size_t detect
  (
    const image::Image<unsigned char> & ima,
    std::unique_ptr<features::Regions> & regions_to_track,
    const size_t min_count,
    const size_t max_count
  ) const override
  {
    // Cast region
    regions_to_track.reset(new RegionT);
    RegionT * regionsCasted = dynamic_cast<RegionT*>(regions_to_track.get());

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
    // Return number of detected features
    return regionsCasted->Features().size();
  }

  bool allocate
  (
    std::unique_ptr<features::Regions> & regions,
    size_t max_elements
  ) override
  {
    regions.reset(new RegionT());
    RegionT * regionsCasted = dynamic_cast<RegionT*>(regions.get());
    regionsCasted->Features().reserve(max_elements);
    regionsCasted->Descriptors().reserve(max_elements);
    return true;
  }

  bool resize
  (
    features::Regions * regions,
    size_t n_elements
  ) override
  {
    RegionT * regionsCasted = dynamic_cast<RegionT *>(regions);
    regionsCasted->Features().resize(n_elements);
    regionsCasted->Descriptors().resize(n_elements);
    return true;
  }

  bool describe
  (
    const image::Image<unsigned char> & ima,
    features::Regions * regions
  ) override
  {
    RegionT * regionsCasted = dynamic_cast<RegionT *>(regions);
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
    return true;
  }
/*
  bool describe
  (
    const image::Image<unsigned char> & ima,
    typename RegionT::FeatureT & new_pt,
    void * desc_raw
  ) override
  {
    RegionT::DescriptorT::bin_type * desc_raw_casted = static_cast<RegionT::DescriptorT::bin_type *>(desc_raw);
    features::PickASDipole(ima, new_pt.x(), new_pt.y(), 10.5f, 0.0f,desc_raw_casted);
    return true;
  }*/

  void getDescriptorRaw
  (
    features::Regions * regions,
    const size_t i,
    void ** desc
  ) override
  {
    //RegionT * regionsCasted = dynamic_cast<RegionT*>(regions);
    *desc =  dynamic_cast<RegionT*>(regions)->Descriptors()[i].data();
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
/*

  size_t InsertRegion
  (
    features::Regions * regions,
    RegionT::FeatureT & pt,
    void * desc
  ) override
  {
    RegionT * regionsCasted = dynamic_cast<RegionT *>(regions);
    RegionT::DescriptorT * descCasted = static_cast<RegionT::DescriptorT*>(desc);
    // Copy features
    regionsCasted->Features().emplace_back(pt);
    regionsCasted->Descriptors().emplace_back(*descCasted);
    return regionsCasted->RegionCount()-1;
  }*/

};

}
}
