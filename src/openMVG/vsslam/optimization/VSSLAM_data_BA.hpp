
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/sfm/sfm_data_BA.hpp>

using namespace openMVG;

namespace openMVG {
namespace VSSLAM {


class VSSLAM_Bundle_Adjustment
{
  public:
    virtual ~VSSLAM_Bundle_Adjustment() = default;

    // Perform a Bundle Adjustment on the SLAM scene (refinement only asked parameters)
    virtual bool Adjust
    (
      // the SfM scene to refine
      VSSLAM_Data & slam_data,
      // tell which parameter needs to be adjusted
      const sfm::Optimize_Options options
    ) =0;
};

}
}
