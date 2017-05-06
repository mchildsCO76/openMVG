// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <stdlib.h>
#include <unistd.h>
#include "openMVG/types.hpp"
#include <openMVG/vsslam/system/Frame.hpp>
#include "openMVG/matching/indMatchDecoratorXY.hpp"

#ifndef SWINE_NOGL
#include "software/VSSLAM/CGlWindow.hpp" // swine
#endif // !SWINE_NOGL

namespace openMVG {
namespace vsslam {

/// Define struct for all parameters

struct VSSLAM_Time_Stats
{
  bool b_enable_time = true;

  double d_feat_detection = 0.0;
  double d_feat_matching = 0.0;
  double d_pose_init = 0.0;

};


}
}
