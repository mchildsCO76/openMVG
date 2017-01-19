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
  virtual ~Abstract_FeatureExtractor(){};

  virtual size_t detect
  (
    const image::Image<unsigned char> & ima,
    std::unique_ptr<features::Regions> & regions_to_track,
    const size_t min_count,
    const size_t max_count
  ) const = 0;

  virtual bool allocate
  (
    std::unique_ptr<features::Regions> & regions,
    size_t max_elements
  ) = 0;
  virtual bool resize
  (
    features::Regions * regions,
    size_t n_elements
  )=0;
  virtual bool describe
  (
    const image::Image<unsigned char> & ima,
    features::Regions * regions
  )=0;

  virtual void getDescriptorRaw
  (
    features::Regions * regions,
    const size_t i,
    void ** desc
  ) =0;
  virtual double SquaredDescriptorDistance
  (
    void * desc_A,
    void * desc_B
  ) =0;

};

} // namespace VSSLAM
} // namespace openMVG
