
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/geometry/Similarity3.hpp>

namespace openMVG  {
namespace VSSLAM  {

using namespace openMVG;
using namespace openMVG::geometry;

  struct MotionModel
  {
    Similarity3 origin_;
    Similarity3 speed_;

    void setSpeed(Similarity3 & pose_t1, Similarity3 & pose_t2)
    {

    }

  };

}
}
