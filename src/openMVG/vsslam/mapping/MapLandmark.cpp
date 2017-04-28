// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <openMVG/vsslam/mapping/MapLandmark.hpp>

namespace openMVG {
namespace vsslam {

  void MapLandmark::addObservation(Frame * frame, const IndexT & feat_id)
  {
    obs_[frame->getFrameId()] = MapObservation(feat_id,frame);
    updateData();
  }

  bool MapLandmark::hasFrameObservation(const IndexT & frame_id)
  {
    if (obs_.find(frame_id) != obs_.end())
      return true;
    return false;
  }

  void MapLandmark::updateData()
  {
    // Update average scale of the descriptor
    feat_mean_scale_ = 0.0;

    // Compute average of scale
    for (auto & o_i : obs_)
    {
      feat_mean_scale_ += o_i.second.frame_ptr->getFeatureScale(o_i.second.feat_id);
    }
    feat_mean_scale_ = feat_mean_scale_/obs_.size();
  }

}
}
