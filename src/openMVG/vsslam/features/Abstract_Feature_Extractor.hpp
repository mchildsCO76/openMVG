// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/types.hpp>
#include <openMVG/vsslam/system/Frame.hpp>
#include <openMVG/vsslam/vsslam_parameters.hpp>

namespace openMVG {
namespace vsslam {

class Abstract_Feature_Extractor
{
protected:
  /// Parameters Object
  std::shared_ptr<VSSLAM_Parameters> params_;

public:
  // Threshold for descriptor matching
  float f_max_desc_dist_high_ = 0.0f;
  float f_max_desc_dist_low_ = 0.0f;

  Abstract_Feature_Extractor(std::shared_ptr<VSSLAM_Parameters> & params)
  {
    params_ = params->share_ptr();
  }

  virtual ~Abstract_Feature_Extractor(){};

  void setParameters(std::shared_ptr<VSSLAM_Parameters> & params)
  {
    params_ = params->share_ptr();
  }

  virtual size_t detect
  (
    const image::Image<unsigned char> & ima,
    Frame * frame,
    const image::Image<unsigned char> * mask = nullptr
  ) const =0;

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

  virtual double getDescriptorDistanceSquared
  (
    void * desc_A,
    void * desc_B
  ) const =0;


};

}
}
