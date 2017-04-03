
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <tuple>
#include <deque>
#include <openMVG/sfm/sfm_data_BA.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include <openMVG/vsslam/Frame.hpp>

using namespace openMVG;

namespace openMVG {
namespace VSSLAM {


class VSSLAM_Bundle_Adjustment
{

public:
  struct BA_options
  {
    bool b_verbose_;
    MAP_POINT_TYPE map_landmark_type_;
    MAP_CAMERA_TYPE map_camera_type_;

    BA_options
    (
      const MAP_CAMERA_TYPE map_camera_type = MAP_CAMERA_TYPE::GLOBAL,
      const MAP_POINT_TYPE map_landmark_type = MAP_POINT_TYPE::GLOBAL_EUCLIDEAN,
      const bool b_verbose = true
    ): b_verbose_(b_verbose), map_camera_type_(map_camera_type), map_landmark_type_(map_landmark_type){};
  };

  virtual ~VSSLAM_Bundle_Adjustment() = default;

  virtual bool OptimizePose
  (
    Frame * frame_i,
    Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_i_idx,
    bool b_use_loss_function = true
  ) =0;

  virtual bool OptimizePose
  (
    Frame * frame_i,
    bool b_use_loss_function = true
  ) =0;


  virtual bool OptimizeLocal
  (
    Frame * frame_i,
    std::vector<std::unique_ptr<MapLandmark> > & vec_triangulated_pts,
    bool b_use_loss_function = true
  ) =0;

  virtual bool addObservationToGlobalSystem(MapLandmark * map_point, MapObservation * map_observation) =0;
  virtual bool addLandmarkToGlobalSysyem(MapLandmark * map_point) =0;
  virtual bool addFrameToGlobalSystem(Frame * frame, bool b_frame_fixed) =0;
  virtual bool optimizeGlobal(MapFrames & map_frames, MapLandmarks & map_landmarks) =0;
};

}
}
