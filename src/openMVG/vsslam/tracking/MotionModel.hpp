// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/types.hpp"
#include <openMVG/vsslam/system/Frame.hpp>
#include <openMVG/vsslam/optimization/Sim3.hpp>

namespace openMVG {
namespace vsslam {


class MotionModel
{
private:
  Vec7 velocity_;  // w | v | s -> sim3
  bool b_valid_ = false;
public:
  bool isValid()
  {
    return b_valid_;
  }
  void setInvalid()
  {
    b_valid_ = false;
  }
  void update(Frame * frame_prev, Frame * frame_cur)
  {
    double time_step = frame_cur->time() - frame_prev->time();
    // TODO: Check how to calculate if we have relative cameras
    Mat4 T_prev_inv = frame_prev->getTransformationMatrixInverse();  // Tp^-1
    Mat4 T_cur = frame_cur->getTransformationMatrix(); // Tc

    // Transformation between two frames  Tp^-1 * Tc
    Mat4 T_v = T_prev_inv*T_cur;

    // Compensate for time difference
    Sim3_log(T_v, velocity_);
    velocity_ = velocity_ / time_step;
    // Set as valid
    b_valid_ = true;
  }

  Mat4 predict(Frame * frame_prev, Frame * frame_cur)
  {
    double time_step = frame_cur->time() - frame_prev->time();
    Mat4 T_prev = frame_prev->getTransformationMatrix(); // Tp
    Vec7 sim3_prev;
    Sim3_log(T_prev, sim3_prev);

    // Compute current pose in sim3
    Vec7 sim3_cur;
    // Add the velocity vector
    sim3_cur = sim3_prev + (velocity_ * time_step);
    sim3_cur(6) = sim3_prev(6);
    // Back to Sim3
    Mat4 T_cur;
    Sim3_exp(sim3_cur,T_cur);

    return T_cur;
  }

};


}
}
