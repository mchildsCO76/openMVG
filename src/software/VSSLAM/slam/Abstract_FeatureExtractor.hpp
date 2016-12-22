// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <software/VSSLAM/slam/Frame.hpp>

namespace openMVG  {
namespace VSSLAM  {

struct Abstract_FeatureExtractor
{
  virtual bool detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count
  ) const = 0;
};

} // namespace VSSLAM
} // namespace openMVG
