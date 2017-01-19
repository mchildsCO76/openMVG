
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FRAME_VSSLAM_HPP
#define FRAME_VSSLAM_HPP

#include <iostream>
#include <memory>

#include <openMVG/numeric/numeric.h>
#include <openMVG/features/features.hpp>
#include <openMVG/geometry/Similarity3.hpp>
#include "openMVG/multiview/projection.hpp"
#include <openMVG/cameras/cameras.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>

namespace openMVG  {
namespace VSSLAM  {

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;

struct MapLandmark;

/// Frame
class Frame : public std::enable_shared_from_this<Frame>
{
public:
  // Basic stats
  size_t frameId_;

  // Camera
  size_t camId_;
  IntrinsicBase * cam_intrinsic_;

  /// Detected features
  std::unique_ptr<features::Regions> regions_;
  std::vector<Vec2> pts_undist_;
  size_t n_feats_ = 0;

  std::vector<MapLandmark *> map_points_; // NULL pointer means no association

  /// Pose
  Similarity3 pose_;
  Mat34 P_;
  double AC_threshold_ = 0.0f;
  // Owner
  size_t ownerFrameId_ = std::numeric_limits<size_t>::max(); // undefined
  Frame * ref_frame_ = nullptr;


  Frame
  (
    const size_t fId,
    const size_t camId,
    IntrinsicBase * cam_intrinsic
  ): frameId_(fId), camId_(camId), cam_intrinsic_(cam_intrinsic)
  {}

  std::shared_ptr<Frame> share_ptr()
  {
    return shared_from_this();
  }
  size_t & getFrameId()
  {
    return frameId_;
  }
  size_t & getCamId()
  {
    return camId_;
  }

  IntrinsicBase * getCameraIntrinsics()
  {
    return cam_intrinsic_;
  }

  size_t & N_feats()
  {
    return n_feats_;
  }

  Vec2 & getFeaturePosition(const size_t & i)
  {
    return pts_undist_[i];
  }

  void updateFeaturesData(bool bCalibratedCamera = true)
  {
    n_feats_ = regions_->RegionCount();
    map_points_ = std::vector<MapLandmark *>(n_feats_,static_cast<VSSLAM::MapLandmark *>(nullptr));

    pts_undist_.resize(n_feats_);

    // If we have distortion and can model it we remove it for faster processing later
    if (bCalibratedCamera && cam_intrinsic_->have_disto())
    {
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for
      #endif
      for ( int i = 0; i < n_feats_; ++i )
      {
        pts_undist_[i] = cam_intrinsic_->remove_disto(regions_->GetRegionPosition(i));
      }
    }
    else
    {
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for
      #endif
      for ( int i = 0; i < n_feats_; ++i )
      {
        pts_undist_[i] = regions_->GetRegionPosition(i);
      }
    }
    pts_undist_.shrink_to_fit();
  }

  // ----------------------
  // -- POSE: c_T_w
  // ----------------------
  void setPose_Sim(const geometry::Similarity3 &sim_pose)
  {
    pose_ = sim_pose;
    updateProjectionMatrix();
  }

  void setPose_Rt(const Mat3 & R, const Vec3 & t)
  {
    pose_ = geometry::Similarity3(geometry::Pose3(R,-R.transpose()*t),1.0);
    updateProjectionMatrix();
  }
  void setPose_Rts(const Mat3 & R, const Vec3 & t, const double s = 1.0)
  {
    pose_ = geometry::Similarity3(geometry::Pose3(R,-R.transpose()*t),s);
    updateProjectionMatrix();
  }

  void setPose_Rcs(const Mat3 & R, const Vec3 & c, const double s = 1.0)
  {
    pose_ = geometry::Similarity3(geometry::Pose3(R,c),s);
    updateProjectionMatrix();
  }

  void setReferenceFrame_Rts
  (
      Frame * new_ref_frame,
      const Mat3 & R,
      const Vec3 & t,
      const double s = 1.0
  )
  {
    pose_ = Similarity3(Pose3(R,-R.transpose()*t),s);
    updateProjectionMatrix();
    ref_frame_ = new_ref_frame;
  }
  void setReferenceFrame
  (
      Frame * new_ref_frame
  )
  {
    ref_frame_ = new_ref_frame;
  }

  void updateProjectionMatrix()
  {
    const Mat3 & cam_K = dynamic_cast<const Pinhole_Intrinsic*>(cam_intrinsic_)->K();
    P_ = cam_K * HStack(pose_.scale_ * pose_.pose_.rotation(),pose_.pose_.translation());
  }


};


} // namespace VSSLAM
} // namespace openMVG


#endif // FRAME_VSSLAM_HPP}

