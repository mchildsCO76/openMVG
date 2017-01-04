
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <software/VSSLAM/slam/Frame.hpp>
#include <software/VSSLAM/slam/Abstract_FeatureExtractor.hpp>


namespace openMVG  {
namespace VSSLAM  {

struct Feat_Extractor_FastDipole : public Abstract_FeatureExtractor
{
  // suggest new feature point for tracking (count point are kept)
  bool detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count
  ) const override
  {
    // If we detect less than that its a fail
    size_t min_detected_features = 100;

    if (count > 0)
    {
      pt_to_track.reserve(count);
      const std::vector<unsigned char> scores = {30, 20, 10, 5};

      features::PointFeatures feats;
      // use a sequence of scores 'in order to deal with lighting change'
      for (const unsigned char fast_score : scores)
      {
        features::FastCornerDetector fastCornerDetector(9, fast_score);
        fastCornerDetector.detect(ima, feats);

        if (feats.size() > count)
        {
          // shuffle to avoid to sample only in one bucket
          std::random_shuffle(feats.begin(), feats.end());
          feats.resize(count); // cut the array to keep only a given count of features

          pt_to_track.swap(feats);
          break;
        }
        else if (fast_score == scores[scores.size()-1] && feats.size() > min_detected_features)
        {
          pt_to_track.swap(feats);
          break;
        }
        feats.clear();
      }
    }
    else
    {
      // if we put count==0 we want all possibilities
      features::FastCornerDetector fastCornerDetector(9, 5);
      fastCornerDetector.detect(ima, pt_to_track);
    }

    if (pt_to_track.size() < min_detected_features)
    {
      return false; // Cannot compute a sufficient number of points for the given image
    }
    return true;
  }

  bool allocate
  (
    std::shared_ptr<Frame> frame,
    size_t max_elements
  ) override
  {
    frame->regions.reset(new features::FAST_Dipole_Regions());
    features::FAST_Dipole_Regions * regionsCasted = dynamic_cast<features::FAST_Dipole_Regions*>(frame->regions.get());
    regionsCasted->Features().reserve(max_elements);
    regionsCasted->Descriptors().reserve(max_elements);
    return true;
  }

  bool describe
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & new_pts,
    std::unique_ptr<features::Regions> &regions
  ) override
  {
    const size_t n_new_pts = new_pts.size();
    // Create region
    regions.reset(new features::FAST_Dipole_Regions);
    features::FAST_Dipole_Regions * regionsCasted = dynamic_cast<features::FAST_Dipole_Regions*>(regions.get());
    // Resize regions to fit new points
    regionsCasted->Features().resize(n_new_pts);
    regionsCasted->Descriptors().resize(n_new_pts);
    // Copy features
    std::copy(new_pts.begin(),new_pts.end(),regionsCasted->Features().begin());

    // Get descriptors of new points
    for (size_t i = 0; i < new_pts.size(); ++i)
    {
      features::PickASDipole(ima, new_pts[i].x(), new_pts[i].y(), 10.5f, 0.0f, regionsCasted->Descriptors()[i].data());
    }
    return true;
  }

  bool insert
  (
    std::shared_ptr<Frame> current_frame,
    std::unique_ptr<features::Regions> &new_regions
  ) override
  {
    features::FAST_Dipole_Regions * regionsCasted = dynamic_cast<features::FAST_Dipole_Regions*>(current_frame->regions.get());
    features::FAST_Dipole_Regions * newRegionsCasted = dynamic_cast<features::FAST_Dipole_Regions*>(new_regions.get());

    // Copy features
    std::copy(newRegionsCasted->Features().begin(),newRegionsCasted->Features().end(),std::back_inserter(regionsCasted->Features()));
    std::copy(newRegionsCasted->Descriptors().begin(),newRegionsCasted->Descriptors().end(),std::back_inserter(regionsCasted->Descriptors()));

    return true;
  }
};

}
}
