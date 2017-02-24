
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <tuple>
#include <deque>
#include <openMVG/sfm/sfm_data_BA.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>

using namespace openMVG;

namespace openMVG {
namespace VSSLAM {


class VSSLAM_Bundle_Adjustment
{
  public:
    virtual ~VSSLAM_Bundle_Adjustment() = default;
/*
    virtual bool OptimizePose
    (
      std::vector<Frame*> * vec_frames,
      Frame * frame_i,
      Hash_Map<MapLandmark *,IndexT> * matches_3D_pts_frame_i_idx,
      std::vector<std::unique_ptr<MapLandmark> > * vec_triangulated_pts
    ) =0;*/
    virtual bool OptimizePose
    (
      Frame * frame_i,
      Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_i_idx
    ) =0;

    virtual bool OptimizeLocal
    (
      Hash_Map<Frame*, size_t> & tmp_frames,
      Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> > & tmp_structure,
      Frame * frame_i,
      std::vector<std::unique_ptr<MapLandmark> > & vec_triangulated_pts
    ) =0;

};

}
}
