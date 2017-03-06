
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <memory>

#include <openMVG/types.hpp>
#include <openMVG/numeric/numeric.h>
#include <openMVG/features/features.hpp>
#include "openMVG/multiview/projection.hpp"
#include <openMVG/cameras/cameras.hpp>
#include "openMVG/matching/matcher_cascade_hashing.hpp"

#include <openMVG/vsslam/Camera.hpp>
#include <openMVG/vsslam/mapping/MapLandmark.hpp>

namespace openMVG  {
namespace VSSLAM  {

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;

class MapLandmark;
struct MapObservation;

/// Frame
class Frame : public std::enable_shared_from_this<Frame>
{
private:
  // Basic stats
  IndexT frame_id_;

  // Camera
  Camera * cam_;

  // Is in the map -> Active
  bool active_;

  // Number of features detected
  IndexT n_feats_ = 0;
  // Detected features
  std::unique_ptr<features::Regions> regions_;

public:

  // Covariance of each detection
  std::vector<Mat> pts_cov_;
  // Distorted/undistorted points - if camera is calibrated (for faster reading)
  std::vector<Vec2> pts_undist_;

  // Associated features with map points
  std::vector<MapLandmark *> map_points_; // NULL pointer means no association

  // Local connectivity map
  Hash_Map<Frame *, size_t> connected_frames_weight;
  std::vector<Frame *> ordered_connected_frames;
  std::vector<size_t> ordered_connected_frames_weights;

  /// Pose
  /// Transformation between: (c - current frame; r - reference frame {e.g W})
  Mat34 P_cr_;
  Mat34 P_rc_;

  Mat4 T_cr_;
  Mat4 T_rc_;

  Mat3 R_rc_;
  Mat3 R_cr_;

  Vec3 O_w_; // Position of the camera in world coordinates

  double AC_reprojection_thresh_ = 0.0f;
  // Owner
  IndexT owner_frame_id_ = UndefinedIndexT; // undefined
  Frame * ref_frame_ = nullptr;


  Frame
  (
    const IndexT fId,
    const IndexT camId,
    Camera * cam
  );

  std::shared_ptr<Frame> share_ptr()
  {
    return shared_from_this();
  }

  const bool & isActive() const
  {
    return active_;
  }

  void setActive()
  {
    active_ = true;
  }
  void setInactive()
  {
    active_ = false;
  }
  const IndexT & getFrameId() const
  {
    return frame_id_;
  }
  const IndexT & getCamId() const
  {
    return cam_->cam_id;
  }

  const bool & getCamCalibrated() const
  {
    return cam_->bCalibrated;
  }

  IntrinsicBase * & getCameraIntrinsics() const
  {
    return cam_->cam_intrinsic_ptr;
  }

  const IndexT & getNumberOfFeatures() const
  {
    return n_feats_;
  }

  features::Regions * getRegions() const
  {
    return regions_.get();
  }

  std::unique_ptr<features::Regions> & getRegionsRaw()
  {
    return regions_;
  }

  Vec2 getFeaturePositionDetected(const IndexT & i) const
  {
    return regions_->GetRegionPosition(i);
  }

  const Vec2 & getFeaturePosition(const IndexT & i) const
  {
    // If camera is calibrated we precalculate the undistorted points
    return pts_undist_[i];
  }

  void updateFeaturesData();
  void undistortPoints();

  bool isPointInFrame(const Vec2 & pt) const;

  void clearMapPoints();

  void clearMapPoint(const IndexT & p_i);

  void setMapPoint(const IndexT p_i, MapLandmark * ml);

  IndexT getNumberMapPoints() const;

  double getSquaredReprojectionError(const Vec3 & pt_frame, const IndexT feat_id) const;


  // ----------------------
  // -- Set reference frame:
  // -- Pose: Transformation from reference (r) to current (c) frame (c_T_r)
  // ----------------------
  void setReferenceFrame (Frame * new_ref_frame);

  void setReferenceFrame_cr_T (Frame * new_ref_frame, const Mat4 &T);

  void setReferenceFrame_cr_Rts
  (
    Frame * new_ref_frame,
    const Mat3 & R,
    const Vec3 & t,
    const double s = 1.0
  );


  void setPose_cr_T
  (
    const Mat4 & this_T_tmp_ref,
    Frame * tmp_ref = nullptr
  );

  void setPose_cr_Rts
  (
    const Mat3 & R,
    const Vec3 & t,
    const double & s,
    Frame * tmp_ref = nullptr
  );

  void getPose_cr_Rts
  (
    Mat3 & R,
    Vec3 & t,
    double & s,
    Frame * tmp_ref = nullptr
  ) const;

  const Vec3 & getCameraCenter() const;



  void updateMatrices();

  Mat34 getProjectionMatrix() const;

  Mat3 getK() const
  {
    return dynamic_cast<const Pinhole_Intrinsic*>(cam_->cam_intrinsic_ptr)->K();;
  }

  Mat3 getK_inv() const
  {
    return dynamic_cast<const Pinhole_Intrinsic*>(cam_->cam_intrinsic_ptr)->Kinv();;
  }

  const Mat4 & getTransformationMatrix_cr() const
  {
    return T_cr_;
  }
  const Mat4 & getTransformationMatrix_rc() const
  {
    return T_rc_;
  }

  const Mat3 & getRotationMatrix_cr() const
  {
    return R_cr_;
  }
  const Mat3 & getRotationMatrix_rc() const
  {
    return R_rc_;
  }



  // ----------------------
  // -- Covisibility stuff
  // ----------------------

  // Compute Frame visibility connections
  void computeFrameVisibilityConnections();

  // Add new Frame connection and update list
  void addFrameVisibilityConnection(Frame * frame, size_t & weight);

  // Update list of best frame visibility connections
  void updateBestFrameVisibilityConnections();

  // Return n-best frames in covisibility graph
  std::vector<Frame *> getBestCovisibilityFrames(const size_t &n);

  // Computes covisibility connecctions on the fly (to prevent screwing graph before putting frame in the system)
  void getFrameVisibilityConnections
  (
    std::vector<Frame *> & temp_ordered_connected_frames,
    const size_t N = std::numeric_limits<size_t>::max()
  ) const;
};


} // namespace VSSLAM
} // namespace openMVG


