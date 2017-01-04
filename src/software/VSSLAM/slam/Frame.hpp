
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FRAME_VSSLAM_HPP
#define FRAME_VSSLAM_HPP

#include <iostream>
#include <memory>

#include <openMVG/features/features.hpp>

namespace openMVG  {
namespace VSSLAM  {

/// Frame
struct Frame:enable_shared_from_this<Frame>
{
  Frame(const uint32_t fId): frameId_(fId) {}

  std::shared_ptr<Frame> share_ptr()
  {
    return shared_from_this();
  }

  uint32_t frameId_;

  /// Detected features
  size_t n_feats;
  std::unique_ptr<features::Regions> regions;
};


} // namespace VSSLAM
} // namespace openMVG


#endif // FRAME_VSSLAM_HPP}
