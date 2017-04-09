
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/vsslam/optimization/sim3.hpp>

namespace openMVG  {
namespace VSSLAM  {

using namespace openMVG;
using namespace openMVG::geometry;

  struct MotionModel
  {
    // velocity - in robotics field w_T_c
    Eigen::Vector7d velocity_;  // w | v | s -> sim3

    Mat4 T_velocity_;

    bool b_valid_ = false;

    bool isValid()
    {
      return b_valid_;
    }
    void setInvalid()
    {
      b_valid_ = false;
    }
    void updateMotionModel(Frame * prev_frame, Frame * cur_frame)
    {
      double time_step = cur_frame->timestamp_ - prev_frame->timestamp_;
      // TODO: Check how to calculate if we have relative cameras
      Mat4 T_prev_inv = prev_frame->getTransformationMatrix();  // Tp^-1
      Mat4 T_cur = cur_frame->getTransformationMatrixInverse(); // Tc

      // Transformation between two frames  Tp^-1 * Tc
      Mat4 T_v = T_prev_inv*T_cur;

      // Compensate for time difference
      Sim3_log(T_v, velocity_);
      velocity_ = velocity_ / time_step;
      // Set as valid
      b_valid_ = true;
    }

    Mat4 predict(Frame * prev_frame, Frame * cur_frame)
    {
      double time_step = cur_frame->timestamp_ - prev_frame->timestamp_;
      Mat4 T_prev = prev_frame->getTransformationMatrixInverse(); // Tp
      Eigen::Vector7d sim3_prev;
      Sim3_log(T_prev, sim3_prev);

      // Compute current pose in sim3
      Eigen::Vector7d sim3_cur;
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
