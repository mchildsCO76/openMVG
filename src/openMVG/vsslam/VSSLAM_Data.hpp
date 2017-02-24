
// Copyright (c) 2017 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/types.hpp"

#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/MapLandmark.hpp>


using namespace openMVG;

namespace openMVG  {
namespace VSSLAM  {


enum MAP_POINT_TYPE
{
  GLOBAL_EUCLIDEAN = 0,  // Absolute world reference frame
  LOCAL_INV_DEPTH = 1,  // Relative reference camera reference frame
};
enum MAP_CAMERA_TYPE
{
  GLOBAL = 0, // Absolute world reference frame
  RELATIVE = 1  // Relative reference camera reference frame
};

/// Define a collection of Pose (indexed by frame::frameId_)
using MapFrames = Hash_Map<IndexT, std::shared_ptr<Frame> >;

/// Define a collection of IntrinsicParameter (indexed by id_intrinsic)
using Intrinsics = Hash_Map<IndexT, cameras::IntrinsicBase *>;


/// Define a collection of Landmarks of the map (3D reconstructed points)
using MapLandmarks = Hash_Map<IndexT,std::unique_ptr<MapLandmark> >;


}
}
