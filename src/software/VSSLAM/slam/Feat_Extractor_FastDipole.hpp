
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

  // suggest new feature point for tracking (count point are kept)
  size_t detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count
  ) const override
  {
    // If we detect less than that its a fail
    size_t min_detected_features = 10;

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
        else if (fast_score == scores[scores.size()-1])
        {
          // if it is the last step we return what we have
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

    // Return number of detected features
    return pt_to_track.size();
  }

  bool allocate
  (
    std::unique_ptr<features::Regions> & regions,
    size_t max_elements
  ) override
  {
    regions.reset(new features::FAST_Dipole_Regions());
    features::FAST_Dipole_Regions * regionsCasted = dynamic_cast<features::FAST_Dipole_Regions*>(regions.get());
    regionsCasted->Features().reserve(max_elements);
    regionsCasted->Descriptors().reserve(max_elements);
    return true;
  }

  bool resize
  (
    std::unique_ptr<features::Regions> & regions,
    size_t n_elements
  ) override
  {
    regions.reset(new features::FAST_Dipole_Regions());
    features::FAST_Dipole_Regions * regionsCasted = dynamic_cast<features::FAST_Dipole_Regions*>(regions.get());
    regionsCasted->Features().resize(n_elements);
    regionsCasted->Descriptors().resize(n_elements);
    return true;
  }

  bool describe
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & new_pts,
    std::unique_ptr<features::Regions> & regions
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


  bool getDescriptorFromFrame
  (
    const std::shared_ptr<Frame> & ref_frame,
    const size_t feature_id,
    openMVG::Vecf * desc
  )
  {
    /*if (feature_id < ref_frame->regions->RegionCount())
    {
      features::FAST_Dipole_Regions * regionsCasted = dynamic_cast<features::FAST_Dipole_Regions*>(ref_frame->regions.get());
      desc = &(regionsCasted->Descriptors()[feature_id]);
      return true;
    }*/
    return false;
  }


  double SquaredDescriptorDistance
  (
      openMVG::Vecf * desc_A,
      openMVG::Vecf * desc_B
  )
  {
    /*
    features::FAST_Dipole_Regions::DescriptorT* d_A = dynamic_cast<features::FAST_Dipole_Regions::DescriptorT*>(desc_A);
    features::FAST_Dipole_Regions::DescriptorT* d_B = dynamic_cast<features::FAST_Dipole_Regions::DescriptorT*>(desc_B);
    openMVG::matching::L2_Vectorized<features::FAST_Dipole_Regions::FeatureT> metric;
    return metric(d_A->data(), d_B->data(), features::FAST_Dipole_Regions::DescriptorT::static_size);*/
    return 50.0;
  }


  bool insert
  (
    std::unique_ptr<features::Regions> & regions,
    features::PointFeature & pt,
    std::unique_ptr<openMVG::Vecf> & desc,
    size_t idx
  ) override
  {
    /*features::FAST_Dipole_Regions * regionsCasted = dynamic_cast<features::FAST_Dipole_Regions*>(current_frame->regions.get());
    features::FAST_Dipole_Regions::DescriptorT descCasted = dynamic_cast<features::FAST_Dipole_Regions::DescriptorT*>(desc.get());
    // Copy features
    regionsCasted->Features()[idx] = pt;
    regionsCasted->Descriptors()[idx] = features::FAST_Dipole_Regions::DescriptorT(*descCasted);*/
    return true;
  }
/*
  size_t getFrameMatching
  (
      std::shared_ptr<Frame> ref_frame,
      const image::Image<unsigned char> & current_ima,
      std::vector<features::PointFeature> & candidate_pts,
      std::vector<bool> & candidate_pts_used,
      const size_t window,
      const float ratio,
      std::unique_ptr<features::Regions> & new_feat_regions,
      Hash_Map<size_t,size_t> & feat_cur_new_matches_ids
  )override
  {
    // Loop through all the features from prev frame
    // and identify the ones in certain window (later checked with descriptors)
    {

      // Extract descriptors that we need
      // vector of pointers (we will use only the ones that are actually used)
      std::vector<std::unique_ptr<features::Descriptor<float, 20> > > candidate_desc(candidate_pts.size());

      std::cout<<"Extracting possible candidate descriptors\n";
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for (int c_i=0; c_i < (int)candidate_pts.size(); ++c_i)
      {
        if (candidate_pts_used[c_i])
        {
          candidate_desc[c_i].reset(new features::Descriptor<float, 20>);
          features::PickASDipole(current_ima, candidate_pts[c_i].x(), candidate_pts[c_i].y(), 10.5f, 0.0f, candidate_desc[c_i]->data());
        }
      }

      size_t valid_desc=0;
      for (int c_i=0; c_i < (int)candidate_pts.size(); ++c_i)
      {
        if(candidate_desc[c_i])
          valid_desc++;
      }
      std::cout<<"Valid desc: "<<valid_desc<<" \n";
      // Compare points

      // Save matches to new_feat_regions


      return 0;
    }
  }
*/

};

}
}
