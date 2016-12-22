// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>

#include <software/VSSLAM/slam/Frame.hpp>

namespace openMVG  {
namespace VSSLAM  {


struct Abstract_Tracker
{
  enum TRACKING_STATUS
  {
    IDLE = 0,
    NOT_INITIALIZED = 0,
    OK = 1,
    LOST = 2
  };

  Abstract_Tracker() = default;

  TRACKING_STATUS tStatus = TRACKING_STATUS::IDLE;
  size_t max_tracked_points = 1500;

  std::shared_ptr<Frame> mPrevFrame;
  std::shared_ptr<Frame> mCurrentFrame;
  std::vector<size_t> prev_tracked_ids;

  // Try to track the pt_to_track features
  virtual bool track
  (
    const image::Image<unsigned char> & ima,
    std::shared_ptr<Frame> current_frame
  ) = 0;

  // suggest new feature point for tracking (count point are kept)
  virtual bool detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count
  ) const = 0;

};

} // namespace VO
} // namespace openMVG
