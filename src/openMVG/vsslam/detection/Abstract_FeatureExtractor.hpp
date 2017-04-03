// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/types.hpp>
#include <openMVG/vsslam/Frame.hpp>

namespace openMVG  {
namespace VSSLAM  {

class Abstract_FeatureExtractor
{
public:
  float max_dist_desc_;
  using RegionT = features::SIFT_Regions;

  virtual ~Abstract_FeatureExtractor(){};

  virtual size_t detect
  (
    const image::Image<unsigned char> & ima,
    Frame * frame,
    const size_t min_count
  ) const = 0;

  virtual bool describe
  (
    const image::Image<unsigned char> & ima,
    const Frame * frame
  ) const =0;

  virtual size_t getDescriptorLength() const =0;
  virtual void getDescriptorRaw
  (
    features::Regions * const regions,
    const IndexT i,
    void ** desc
  ) const =0;
  virtual double SquaredDescriptorDistance
  (
    void * desc_A,
    void * desc_B
  ) const =0;

};

} // namespace VSSLAM
} // namespace openMVG
