
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <tuple>
#include <deque>
#include <openMVG/sfm/sfm_data_BA.hpp>

using namespace openMVG;

namespace openMVG {
namespace VSSLAM {


class VSSLAM_Bundle_Adjustment
{
  public:
    virtual ~VSSLAM_Bundle_Adjustment() = default;
/*
    // Perform a Bundle Adjustment on the SLAM scene (refinement only asked parameters)
    virtual bool Adjust
    (
      // the SfM scene to refine
      VSSLAM_Data & slam_data,
      // tell which parameter needs to be adjusted
      const sfm::Optimize_Options options,
      const bool first_pose_fixed
    ) =0;
*/

    virtual bool OptimizePose
    (
      std::vector<Frame*> & vec_frames,
      Hash_Map<MapLandmark *,size_t> * matches_3D_ptr_cur_idx,
      std::vector<std::pair<Vec3, std::deque<std::pair<Frame*,size_t> > > > * vec_triangulated_pts
    ) =0;

};

}
}
