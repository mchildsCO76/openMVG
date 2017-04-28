// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/types.hpp"
#include "openMVG/cameras/cameras.hpp"

namespace openMVG {
namespace vsslam {

class Frame;
class MapLandmark;

/// Define the frame type representation
enum MAP_FRAME_TYPE
{
  GLOBAL = 0, // Absolute world reference frame
  RELATIVE = 1  // Relative reference camera reference frame
};

/// Define the landmark type representation
enum MAP_LANDMARK_TYPE
{
  GLOBAL_EUCLIDEAN = 0,  // Absolute XYZ
  LOCAL_INV_DEPTH = 1  // Relative with respect to reference frame
};

enum MAP_OPTIMIZATION_TYPE
{
  CERES = 0,
  SLAMPP = 1
};

enum TRACKING_STATUS
{
  NOT_INIT = 0,
  INIT = 1,
  OK = 2,
  LOST = 3,
  IDLE = 4
};

/// Define collection of frames
using MapFrames = Hash_Map<IndexT, std::shared_ptr<Frame> >;
/// Define a collection of landmarks of the map (3D reconstructed points)
using MapLandmarks = Hash_Map<IndexT,std::unique_ptr<MapLandmark> >;
// Define a collection of new landmarks
using NewMapLandmarks = std::vector<std::unique_ptr<MapLandmark> >;

/// Define a collection of IntrinsicParameter (indexed by id_intrinsic)
using Intrinsics = Hash_Map<IndexT, cameras::IntrinsicBase *>;


}
}
