
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching_image_collection/Matcher_Regions_SiftGPU.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matcher_kdtree_flann.hpp"
#include "openMVG/matching/matcher_cascade_hashing.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/matching_image_collection/Matcher.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

namespace openMVG {
namespace matching_image_collection {

using namespace openMVG::matching;
using namespace openMVG::features;

Matcher_Regions_SiftGPU::Matcher_Regions_SiftGPU(
  float distRatio, size_t max_features)
  :Matcher(), _f_dist_ratio(distRatio), _max_features(max_features)
{
	matcher->SetMaxSift(max_features);
}

void Matcher_Regions_SiftGPU::Match(
  const sfm::SfM_Data & sfm_data,
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  const Pair_Set & pairs,
  PairWiseMatches & map_PutativesMatches)const // the pairwise photometric corresponding points
{
	if (regions_provider->regions_per_view.empty())
		return;

	// initialize SiftGPU matcher
	if(matcher->VerifyContextGL() == 0) return;

  C_Progress_display my_progress_bar( pairs.size() );

  // Sort pairs according the first index to minimize the MatcherT build operations
  typedef std::map<size_t, std::vector<size_t> > Map_vectorT;
  Map_vectorT map_Pairs;
  for (Pair_Set::const_iterator iter = pairs.begin(); iter != pairs.end(); ++iter)
  {
    map_Pairs[iter->first].push_back(iter->second);
  }


  // Perform matching between all the pairs
  for (Map_vectorT::const_iterator iter = map_Pairs.begin();
    iter != map_Pairs.end(); ++iter)
  {
    const size_t I = iter->first;
    const std::vector<size_t> & indexToCompare = iter->second;

    const features::SIFT_Regions & regionsI = *(dynamic_cast<SIFT_Regions*>(regions_provider->regions_per_view.at(I).get()));
    if (regionsI.RegionCount() == 0)
    {
      my_progress_bar += indexToCompare.size();
      continue;
    }
    // Initialize the matching interface
    const unsigned char *desc_I = (const unsigned char *)regionsI.DescriptorRawData();
    size_t regions_I_count = regionsI.RegionCount();

    matcher->SetDescriptors(0,regions_I_count,desc_I);

    for (int j = 0; j < (int)indexToCompare.size(); ++j)
    {
      const size_t J = indexToCompare[j];

      const features::SIFT_Regions &regionsJ = *(dynamic_cast<SIFT_Regions*>(regions_provider->regions_per_view.at(J).get()));
      if (regionsJ.RegionCount() == 0
          || regionsI.Type_id() != regionsJ.Type_id())
      {
        ++my_progress_bar;
        continue;
      }

      // Initialize the matching interface
      const unsigned char *desc_J = (const unsigned char *)regionsJ.DescriptorRawData();
      size_t regions_J_count = regionsJ.RegionCount();

      matcher->SetDescriptors(1,regions_J_count,desc_J);
      //matcher->SetDescriptors(1,descriptorsJ.size(),&descJ[0]);

      int (*match_buf)[2] = new int[regions_I_count][2];

      int num_match = matcher->GetSiftMatch(regions_I_count, match_buf,std::numeric_limits<float>::max(),_f_dist_ratio );
      //std::cout << num_match << " sift matches were found;\n";

      IndMatches vec_putatives_matches;
      vec_putatives_matches.reserve(num_match);

      for(int i  = 0; i < num_match; ++i)
      {

    	  vec_putatives_matches.emplace_back(match_buf[i][0],match_buf[i][1]);
      }

      {
        ++my_progress_bar;
        if (!vec_putatives_matches.empty())
        {
          map_PutativesMatches.insert( make_pair( make_pair(I,J), std::move(vec_putatives_matches) ));
        }
      }
    }
  }
}

} // namespace openMVG
} // namespace matching_image_collection
