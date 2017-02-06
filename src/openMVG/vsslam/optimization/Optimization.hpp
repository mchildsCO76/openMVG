
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace openMVG  {
namespace VSSLAM  {

  static bool createInitialPose_VSSLAM_Data
  (
    VSSLAM_Data & init_map,
    std::shared_ptr<Frame> & frame_1,
    std::shared_ptr<Frame> & frame_2,
    Hash_Map<size_t,size_t> & putative_matches_1_2_idx,
    const std::vector<size_t> & inliers
  )
  {
    const size_t frame_id_1 = frame_1->frameId_;
    const size_t frame_id_2 = frame_2->frameId_;
    const Camera * cam_1 = frame_1->cam_;
    const Camera * cam_2 = frame_2->cam_;
    const IntrinsicBase * cam_intrinsic_1 = cam_1->cam_intrinsic_ptr;
    const IntrinsicBase * cam_intrinsic_2 = cam_2->cam_intrinsic_ptr;

    // -------------------
    // -- Construct initial map (only serves for BA)
    // -------------------

    // Add initial keyframes
    init_map.keyframes[frame_id_1] = frame_1->share_ptr();
    init_map.keyframes[frame_id_2] = frame_2->share_ptr();
    init_map.cam_intrinsics[cam_1->cam_id] = cam_1->cam_intrinsic_ptr;
    init_map.cam_intrinsics[cam_2->cam_id] = cam_2->cam_intrinsic_ptr;

    // Get Projection matrices of cameras
    Mat34 P1,P2;

    switch (init_map.mapCameraType)
    {
      case VSSLAM::MAP_CAMERA_TYPE::ABSOLUTE:
        P1 = frame_1->getProjectionMatrix();
        P2 = frame_2->getProjectionMatrix();
      break;
      case VSSLAM::MAP_CAMERA_TYPE::RELATIVE:
      /*{
        // If cameras are relative we find tranformation
        Mat4 T_cam1_2;
        if (!getRelativeCameraTransformation(frame_1.get(), frame_2.get(),T_cam1_2))
        {
          std::cout<<"Cant find relative poses!!\n Something is seriously wrong!\n";
          return false;
        }
        P1 = frame_1->getK()*T_cam1_2.block(0,0,3,4);
        // Set second camera as base
        P2 = HStack(Mat3::Identity(), Vec3::Zero());
      }*/
        return false;
      break;
    }

    // ------------------------
    // Add points to the initial map
    // ------------------------
    for ( size_t k : inliers)
    {
      Hash_Map<size_t,size_t>::const_iterator m_iter = putative_matches_1_2_idx.begin();
      std::advance(m_iter,k);
      const size_t feat_id_1 = m_iter->first;
      const size_t feat_id_2 = m_iter->second;

      // 3D point
      MapLandmark l;
      l.id_ = init_map.getNextFreeLandmarkId();

      // Position of detected features
      //   If we calibrated camera we just use undistorted
      //   If not calibrated we remove the distortion as best as we know it (current estimation)
      const Vec2
        & x1_ = cam_1->bCalibrated ? frame_1->getFeaturePositionUndistorted(feat_id_1) : cam_intrinsic_1->remove_disto(frame_1->getFeaturePositionUndistorted(feat_id_1)),
        & x2_ = cam_2->bCalibrated ? frame_2->getFeaturePositionUndistorted(feat_id_2) : cam_intrinsic_2->remove_disto(frame_2->getFeaturePositionUndistorted(feat_id_2));

      // Triangulate results
      TriangulateDLT(P1, x1_, P2, x2_, &l.X_);

      // Add observations
      MapObservation m1(feat_id_1,frame_1.get());
      MapObservation m2(feat_id_2,frame_2.get());
      // we input distorted/undistorted coordinates based on the camera calibration
      m1.pt_ptr = &(frame_1->getFeaturePositionUndistorted(feat_id_1));
      m2.pt_ptr = &(frame_2->getFeaturePositionUndistorted(feat_id_2));

      l.obs_[frame_id_1] = m1;
      l.obs_[frame_id_2] = m2;

      // Add points to initial map
      init_map.structure[l.id_] = l;

      /*
      if (init_map.mapCameraType == VSSLAM::MAP_CAMERA_TYPE::RELATIVE)
      {
        l.ref_frame_ = frame_2.get();
      }*/
    }
  }
}
}
