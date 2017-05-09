// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include <openMVG/types.hpp>
#include <openMVG/features/regions.hpp>
#include <openMVG/vsslam/system/Camera.hpp>
#include <openMVG/vsslam/optimization/PoseEstimator.hpp>
#include <openMVG/vsslam/optimization/Sim3.hpp>
#include <openMVG/vsslam/mapping/MapLandmark.hpp>

namespace openMVG {
namespace vsslam {


class Frame: public std::enable_shared_from_this<Frame>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW // swine - for some reason does not work with std::shared_ptr :(

private:
  IndexT id_;
  double d_timestamp_;
  bool b_active_; // Flag if frame is in the global map
  Camera * ptr_cam_;  // Pointer to camera object

	typedef Eigen::Matrix<double, 4, 4, Eigen::DontAlign>  Mat4U;
	typedef Eigen::Matrix<double, 3, 3, Eigen::DontAlign>  Mat3U;
	typedef Eigen::Matrix<double, 3, 4, Eigen::DontAlign>  Mat34U; // swine - a quick hack but works

  // Pose
  Mat4U T_;
  Mat4U T_inv_;
  Mat3U R_;
  Mat3U R_inv_;
  Mat34U P_; // Projection matrix (world -> camera) .. associated with T_inv_
  Vec3 origin_;
  Frame * frame_reference_ = nullptr; // owner frame of current frame (nullptr -> global)

  // Scene
  float f_scene_median_;

  // Features
  size_t n_feats_ = 0;  // Number of detected features
  std::unique_ptr<features::Regions> regions_;

  std::vector<Vec2> vec_pts_undist_;  // Vector of undistorted positions of keypoints
  std::vector<Eigen::Matrix<double, 2, 2, Eigen::DontAlign> > vec_pts_sqrt_inf_mat_;  // Vector of information matrix for each keypoint
  std::vector<float> vec_pts_scale_;

  // Landmarks
  std::vector<MapLandmark *> map_landmark_;

public:
  Frame
  (
    const IndexT & id_frame,
    const double & timestamp,
    Camera * ptr_cam
  );

  std::shared_ptr<Frame> share_ptr()
  {
    return shared_from_this();
  }

  const IndexT & getFrameId() const
  {
    return id_;
  }
  const IndexT & getCamId() const
  {
    return ptr_cam_->getCamId();
  }

  const bool & isActive() const
  {
    return b_active_;
  }

  void setActive()
  {
    b_active_ = true;
  }

  void setInactive()
  {
    b_active_ = false;
  }

  const double & time() const
  {
    return d_timestamp_;
  }

  const bool & isCamCalibrated() const
  {
    return ptr_cam_->isCalibrated();
  }

  IntrinsicBase * & getCameraIntrinsics() const
  {
    return ptr_cam_->ptr_intrinsic_valid_;
  }

  const Mat3 & getK() const
  {
    return dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->K();
  }

  const Mat3 & getKinv() const
  {
    return dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->Kinv();
  }


  // ------------------------------------------------
  //  Scene
  // ------------------------------------------------

  const float & getSceneMedian()
  {
    return f_scene_median_;
  }
  // ------------------------------------------------
  //  Features
  // ------------------------------------------------

  features::Regions * getRegions() const
  {
    return regions_.get();
  }

  std::unique_ptr<features::Regions> & getRegionsPtr()
  {
    return regions_;
  }

  void updateFeaturesData();

  void undistortPoints();

  const size_t & getNumberOfFeatures() const
  {
    return n_feats_;
  }
  std::vector<Vec2> & getFeaturePointVector()
  {
    return vec_pts_undist_;
  }
  const std::vector<Vec2> & getFeaturePointVector() const
  {
    return vec_pts_undist_;
  }
  std::vector<float> & getFeatureScaleVector()
  {
    return vec_pts_scale_;
  }
  const std::vector<float> & getFeatureScaleVector() const
  {
    return vec_pts_scale_;
  }
  std::vector<Eigen::Matrix<double, 2, 2, Eigen::DontAlign> > & getFeatureSqrtInfMatrixVector()
  {
    return vec_pts_sqrt_inf_mat_;
  }
  const std::vector<Eigen::Matrix<double, 2, 2, Eigen::DontAlign> > & getFeatureSqrtInfMatrixVector() const
  {
    return vec_pts_sqrt_inf_mat_;
  }

  const Vec2 & getFeaturePosition(const size_t & i) const
  {
    return vec_pts_undist_[i];
  }

  Vec2 & getFeaturePosition(const size_t & i)
  {
    return vec_pts_undist_[i];
  }

  Vec3 getFeaturePositionHomogeneous(const size_t & feat_i)
  {
    return (ptr_cam_->isCalibrated() ? getFeaturePosition(feat_i) : getCameraIntrinsics()->remove_disto(getFeaturePosition(feat_i))).homogeneous();
  }

  bool getProjectedPoint(const Vec3 & pt_3D, const Frame * pt_3D_reference, Vec2 & pt_2D_frame) const;

  bool getProjectedPoint(const MapLandmark * map_landmark, Vec2 & pt_2D_frame) const;

  Vec3 getRayToMapLandmark(const MapLandmark * map_landmark);
  Vec3 getRayToPoint(const Vec3 & pt_3D, const Frame * pt_3D_reference);
  Vec3 getRayToPoint(const IndexT & feat_id);

  const float & getFeatureScale(const size_t & i) const
  {
    return vec_pts_scale_[i];
  }

  float & getFeatureScale(const size_t & i)
  {
    return vec_pts_scale_[i];
  }

  const Eigen::Matrix<double, 2, 2, Eigen::DontAlign> & getFeatureSqrtInfMatrix(const size_t & i) const
  {
    return vec_pts_sqrt_inf_mat_[i];
  }

  Eigen::Matrix<double, 2, 2, Eigen::DontAlign> & getFeatureSqrtInfMatrix(const size_t & i)
  {
    return vec_pts_sqrt_inf_mat_[i];
  }

  // ------------------------------------------------
  //  Feature association
  // ------------------------------------------------
  bool isPointInFrame(const Vec2 & pt) const;
  bool checkFeatureAssociation(const Vec3 & pt3D_frame, const IndexT feat_id, double chi2_px2_thresh) const;
  bool checkFeatureAssociation(const Vec2 & pt_2D_frame, const IndexT feat_id, double thresh) const;
  bool checkLandmarkPosition(const Vec3 & pt3D_frame) const;


  // ------------------------------------------------
  //  Map Landmarks
  // ------------------------------------------------
  std::vector<MapLandmark *>  & getLandmarks()
  {
    return map_landmark_;
  }
  const std::vector<MapLandmark *>  & getLandmarks() const
  {
    return map_landmark_;
  }
  MapLandmark * & getLandmark(IndexT & feat_id)
  {
    return map_landmark_[feat_id];
  }
  void setLandmark(IndexT & feat_id, MapLandmark * map_landmark);
  void removeLandmark(IndexT & feat_id);

  const float getMapLandmarkDepth(MapLandmark * map_landmark) const;

  size_t getNumberOfMapPoints(bool b_only_global = true) const;

  void computeSceneStatistics();


  // ------------------------------------------------
  //  Pose
  // ------------------------------------------------
  Frame * getReferenceFrame() const
  {
    return frame_reference_;
  }
  void setPose_T(const Mat4 & T, Frame * frame_reference);
  void setPose_sim3
  (
    Eigen::Matrix<double, 7, 1> & vec_state,
    Frame * frame_reference
  );
  void setPose_sRt_Inverse(const Mat3 & R, const Vec3 & t, const double & s,Frame * frame_reference);

  void setPose_T_withReferenceFrame (const Mat4 & T, Frame * frame_new_reference);

  void updatePoseData(bool b_inverse_updated = false);



  // Get camera center expressed in world coordinates (regardless of representation of camera)
  const Vec3 getCameraCenter() const;


  const Mat4U & getTransformationMatrix() const
  {
    return T_;
  }
  const Mat4U & getTransformationMatrixInverse() const
  {
    return T_inv_;
  }

  const Mat3U & getRotationMatrix() const
  {
    return R_;
  }
  const Mat3U & getRotationMatrixInverse() const
  {
    return R_inv_;
  }

  const double getPoseScale() const
  {
    return T_.block(0,0,3,1).norm();
  }

  void getPose_T
  (
    Mat4 & T,
    Frame * frame_reference
  ) const;

  void getPose_StateVector
  (
    Eigen::Matrix<double, 12, 1> & vec_state,
    Frame * frame_reference
  ) const;
  void getPose_sRt_Inverse
  (
      Mat3 & R,
      Vec3 & t,
      double & s,
      Frame * frame_reference
  ) const;

  Mat34 getCameraProjectionMatrix(Frame * frame_reference);


  // ------------------------------------------------
  //  Visibility connections - local maps
  // ------------------------------------------------
  void getFrameVisibilityConnections(std::vector<Frame *> & frames_connected_ordered, const size_t n_best) const;

};

}
}
