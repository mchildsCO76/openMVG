
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
#include <openMVG/vsslam/optimization/sim3.hpp>

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

  // Distorted/undistorted points - if camera is calibrated (for faster reading)
  std::vector<Vec2> pts_undist_;
  // Information matrix of each detection
  std::vector<Eigen::Matrix2d> pts_information_mat_;
  // Scale of features
  std::vector<float> pts_scale_;

  // Associated features with map points
  std::vector<MapLandmark *> map_points_; // NULL pointer means no association

  // Local connectivity map
  Hash_Map<Frame *, size_t> connected_frames_weight;
  std::vector<Frame *> ordered_connected_frames;
  std::vector<size_t> ordered_connected_frames_weights;

  /// c - current frame; r - reference frame {e.g W}
  /// transformation from world to camera X_c = T_cr_ * X_w (computer vision)
  /// transformation from camera to world X_r = T_rc_ * X_c (robotics)

  // Projection matrix
  Mat34 P_cr_;
  Mat34 P_rc_; // Do we actually need it?

  // Transformation matrix
  Mat4 T_cr_; // from world to camera (computer vision)
  Mat4 T_rc_; // robotics
  // Rotation matrix
  Mat3 R_cr_;
  Mat3 R_rc_;
  // Position of camera in world coordinates
  Vec3 O_w_;


  // Scene detainls
  float f_scene_median = 0.0;
  float f_scene_095 = std::numeric_limits<float>::max();
  double timestamp_;


  double reproj_thresh_sq_ = 4.0f;
  // Owner
  IndexT owner_frame_id_ = UndefinedIndexT; // undefined
  Frame * ref_frame_ = nullptr; // owner frame of current frame (nullptr -> global)


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

  const float & getFeatureScale(const IndexT & i) const
  {
    // If camera is calibrated we precalculate the undistorted points
    return pts_scale_[i];
  }

  const float getFeatureSigmaSq(const IndexT & i) const
  {
    // If camera is calibrated we precalculate the undistorted points
    return 16.0f; //4.0f*4.0f
  }

  const Eigen::Matrix2d & getFeatureInformationMatrix(const IndexT & i) const
  {
    // Return information matrix of the current feature
    return pts_information_mat_[i];
  }

  void updateFeaturesData();
  void undistortPoints();

  bool isPointInFrame(const Vec2 & pt) const;

  void clearMapPoints();

  void clearMapPoint(const IndexT & p_i);

  void setMapPoint(const IndexT p_i, MapLandmark * ml);

  size_t getNumberMapPoints() const;
  size_t getNumberGlobalMapPoints() const;

  double getSquaredReprojectionError(const Vec3 & pt_frame, const IndexT feat_id) const;

  bool checkFeatureAssociation(const Vec3 & pt3D_frame, const IndexT feat_id, double chi2_px2_thresh) const;

  void computeSceneStatistics();

  bool checkLandmarkPosition(const Vec3 & map_landmark_3D);
  /// c - current frame; r - reference frame {e.g W}
  /// transformation from world to camera X_c = T_cr_ * X_w (computer vision)
  /// transformation from camera to world X_r = T_rc_ * X_c (robotics)

  // Update all pose related matrices
  // Everything is updated from T_cr_
  void updatePoseMatrices();

  // Get pose of camera in desired (tmp_ref) reference frame
  // aka get transformation from desired reference (tmp_ref) frame to this frame
  // (e.g. from world to camera reference frame)
  void getPose_Rts
  (
    Mat3 & R,
    Vec3 & t,
    double & s,
    Frame * tmp_ref // Pose [sR t] expressed in tmp_ref reference frame
  ) const;

  // Get transformation from this frame to desired reference (tmp_ref) frame
  // (e.g. from camera to world reference frame)
  void getPoseInverse_Rts
  (
    Mat3 & R,
    Vec3 & t,
    double & s,
    Frame * tmp_ref // Pose [sR t] expressed in tmp_ref reference frame
  ) const;


  // Get transformation from this frame to desired reference (tmp_ref) frame expressed in sim3
  // (e.g. from camera to world reference frame)
  void getPoseInverse_sim3
  (
    Eigen::Matrix<double, 7, 1> & v_state,
    Frame * tmp_ref
  ) const;

  // Get transformation from this frame to desired reference (tmp_ref) frame expressed with state vector
  // (e.g. from camera to world reference frame)
  void getPoseInverse_StateVector
  (
    Eigen::Matrix<double, 12, 1> & v_state,
    Frame * tmp_ref // Pose [sR t] expressed in tmp_ref reference frame
  ) const;

  // Get camera center expressed in world coordinates (regardless of representation of camera)
  const Vec3 & getCameraCenter() const;


  // Set pose of the frame with transformation T expressed
  // as transformation from tmp_reference to this frame (camera) : c_T_tmp_ref
  // If reference frame of the given transformation (tmp_ref) is not equal to reference frame of the frame (ref_frame)
  // we need transformation from reference of the frame to reference of given transformation: tmp_ref_T_ref
  // Then we can compute pose: c_T_ref = c_T_tmp_ref * tmp_ref_T_ref
  void setPose_T
  (
    const Mat4 & c_T_tmp_ref,
    Frame * tmp_ref
  );

  // Set pose of the frame with [sR t] transformation expressed
  // as [sR t] from tmp_reference to this frame (camera) : c_T_tmp_ref=[sR t]
  // If reference frame of the given transformation (tmp_ref) is not equal to reference frame of the frame (ref_frame)
  // we need transformation from reference of the frame to reference of given transformation: tmp_ref_T_ref
  // Then we can compute pose: c_T_ref = c_T_tmp_ref * tmp_ref_T_ref
  void setPose_Rts
  (
    const Mat3 & R,
    const Vec3 & t,
    const double & s,
    Frame * tmp_ref
  );


  // Set pose of the frame with [sR t] transformation expressed
  // as sim3 vector [upsilon, omega, sigma] from this frame (camera) to reference frame : tmp_ref_T_c
  // TODO: If reference frame of the given transformation (tmp_ref) is not equal to reference frame of the frame (ref_frame)
  void setPoseInverse_sim3
  (
    const Eigen::Matrix<double, 7, 1> & v_state,
    Frame * tmp_ref
  );


  // Set frame as the reference frame
  void setReferenceFrame (Frame * new_ref_frame);

  void setReferenceFrameAndPose_T (Frame * new_ref_frame, const Mat4 &T);

  void setReferenceFrameAndPose_Rts
  (
    Frame * new_ref_frame,
    const Mat3 & R,
    const Vec3 & t,
    const double s = 1.0
  );




  Mat34 getProjectionMatrix() const;

  Mat3 getK() const
  {
    return dynamic_cast<const Pinhole_Intrinsic*>(cam_->cam_intrinsic_ptr)->K();;
  }

  Mat3 getK_inv() const
  {
    return dynamic_cast<const Pinhole_Intrinsic*>(cam_->cam_intrinsic_ptr)->Kinv();;
  }

  const Mat4 & getTransformationMatrix() const
  {
    return T_cr_;
  }
  const Mat4 & getTransformationMatrixInverse() const
  {
    return T_rc_;
  }

  const Mat3 & getRotationMatrix() const
  {
    return R_cr_;
  }
  const Mat3 & getRotationMatrixInverse() const
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


