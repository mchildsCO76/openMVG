
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
#include <openMVG/vsslam/Camera.hpp>

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
  Camera * cam_;

  /// Detected features
  std::unique_ptr<features::Regions> regions_;
  std::vector<Vec2> pts_undist_;
  size_t n_feats_ = 0;

  std::vector<MapLandmark *> map_points_; // NULL pointer means no association

  Hash_Map<Frame *, size_t> connected_frames_weight;
  std::vector<Frame *> ordered_connected_frames;
  std::vector<size_t> ordered_connected_frames_weights;


  /// Pose
  Similarity3 pose_;
  Mat34 P_;
  Mat4 T_;
  double AC_reprojection_thresh_ = 0.0f;
  // Owner
  size_t ownerFrameId_ = std::numeric_limits<size_t>::max(); // undefined
  Frame * ref_frame_ = nullptr;


  Frame
  (
    const size_t fId,
    const size_t camId,
    Camera * cam
  ): frameId_(fId), camId_(camId), cam_(cam)
  {
    T_ = Mat4::Identity();
  }

  std::shared_ptr<Frame> share_ptr()
  {
    return shared_from_this();
  }
  const size_t & getFrameId()
  {
    return frameId_;
  }
  size_t & getCamId()
  {
    return camId_;
  }

  IntrinsicBase * getCameraIntrinsics()
  {
    return cam_->cam_intrinsic_ptr;
  }

  size_t & N_feats()
  {
    return n_feats_;
  }

  Vec2 getFeaturePosition(const size_t & i)
  {
    return regions_->GetRegionPosition(i);
  }

  Vec2 & getFeaturePositionUndistorted(const size_t & i)
  {
    // If camera is calibrated we precalculate the undistorted points
    return pts_undist_[i];
  }
  void updateFeaturesData()
  {
    n_feats_ = regions_->RegionCount();
    map_points_ = std::vector<MapLandmark *>(n_feats_,static_cast<VSSLAM::MapLandmark *>(nullptr));
    // If camera is calibrated we remove distortion and save points
//    if (cam_->bCalibrated)
    //    {
      pts_undist_.resize(n_feats_);
      undistortPoints();
      pts_undist_.shrink_to_fit();
      // }

  }

  void clearMapPoints()
  {
    map_points_ = std::vector<MapLandmark *>(n_feats_,static_cast<VSSLAM::MapLandmark *>(nullptr));
  }

  void undistortPoints()
  {
    // If we have distortion and can model it we remove it for faster processing later
    if (cam_->cam_intrinsic_->have_disto() && cam_->bCalibrated)
    {
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for
      #endif
      for ( int i = 0; i < n_feats_; ++i )
      {
        pts_undist_[i] = cam_->cam_intrinsic_->remove_disto(regions_->GetRegionPosition(i));
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
  }

  // ----------------------
  // -- POSE: c_T_w
  // ----------------------
  void setPose_Sim(const geometry::Similarity3 &sim_pose)
  {
    pose_ = sim_pose;
    std::cout<<"Set pose SIM: \nR: "<<sim_pose.rotation()<<"\n t: "<<sim_pose.translation()<<"\n s: "<<sim_pose.scale_<<"\n";
    updateMatrices();
  }

  void setPose_Rts(const Mat3 & R, const Vec3 & t, const double s = 1.0)
  {
    pose_ = geometry::Similarity3(geometry::Pose3(R,-R.transpose()*t),s);
    updateMatrices();
  }

  void setPose_Rcs(const Mat3 & R, const Vec3 & c, const double s = 1.0)
  {
    pose_ = geometry::Similarity3(geometry::Pose3(R,c),s);
    updateMatrices();
  }

  void setPose_T(const Mat4 &T)
  {
    const double scale = T.block(0,0,3,1).norm();
    const Mat3 R = T.block(0,0,3,3)/scale;
    const Vec3 t = T.block(0,3,3,1)/scale;
    pose_ = Similarity3(Pose3(R,-R.transpose()*t),scale);
    updateMatrices();
  }

  void setReferenceFrame_Sim
  (
      Frame * new_ref_frame,
      const geometry::Similarity3 &sim_pose
  )
  {
    ref_frame_ = new_ref_frame;
    pose_ = sim_pose;
    std::cout<<"Set pose SIM: \nR: "<<sim_pose.rotation()<<"\n t: "<<sim_pose.translation()<<"\n s: "<<sim_pose.scale_<<"\n";
    updateMatrices();

  }
  void setReferenceFrame_Rts
  (
      Frame * new_ref_frame,
      const Mat3 & R,
      const Vec3 & t,
      const double s = 1.0
  )
  {
    ref_frame_ = new_ref_frame;
    pose_ = Similarity3(Pose3(R,-R.transpose()*t),s);
    updateMatrices();
  }


  void setReferenceFrame
  (
      Frame * new_ref_frame
  )
  {
    ref_frame_ = new_ref_frame;
  }

  void updateMatrices()
  {
    const Mat3 & cam_K = dynamic_cast<const Pinhole_Intrinsic*>(cam_->cam_intrinsic_ptr)->K();
    P_ = cam_K * HStack(pose_.scale_ * pose_.pose_.rotation(),pose_.scale_ * pose_.pose_.translation());
    T_.block(0,0,3,4) = HStack(pose_.scale_ * pose_.pose_.rotation(),pose_.scale_ * pose_.pose_.translation());
  }

  Mat34 getProjectionMatrix()
  {
    return P_;
  }
  Mat3 getK()
  {
    return dynamic_cast<const Pinhole_Intrinsic*>(cam_->cam_intrinsic_ptr)->K();;
  }

  Mat4 getTransformationMatrix() const
  {
    return T_;
  }
  Mat4 getInverseTransformationMatrix() const
  {
    return T_.inverse();
  }

  Mat3 getInverseRotationMatrix() const
  {
    return 1/pose_.scale_ * pose_.rotation().transpose();
  }

  // ----------------------
  // -- Covisibility stuff
  // ----------------------

  void updateFrameVisibilityConnections()
  {
    std::cout<<"Update Frame covisibility: "<<frameId_<<"\n";
    // Loop through all matched points of the frame and vote for each frame that sees any of the points
    connected_frames_weight.clear();
    for (MapLandmark * map_point: map_points_)
    {
      if (!map_point)
        continue;

      LandmarkObservations & map_obs = map_point->obs_;
      for (auto obs: map_obs)
      {
        if (obs.first == frameId_)
          continue;
        connected_frames_weight[obs.second.frame_ptr]++;
      }
    }
    if (connected_frames_weight.empty())
      return;

    // Add all connections to frames that have at least min_points in common
    int min_common_pts = 15;
    Frame * bestFrame = nullptr;
    size_t max_common_pts = 0;

    std::vector<std::pair<size_t,Frame*> > vec_frames_pts;
    vec_frames_pts.reserve(connected_frames_weight.size());

    for (auto m_iter : connected_frames_weight)
    {
      size_t & n_common_pts = m_iter.second;
      if (n_common_pts > max_common_pts)
      {
        max_common_pts = n_common_pts;
        bestFrame = m_iter.first;
      }

      if (n_common_pts > min_common_pts)
      {
        vec_frames_pts.emplace_back(std::make_pair(n_common_pts, m_iter.first));
        // Add connection to new frame in the destination frame
        m_iter.first->addFrameVisibilityConnection(this, n_common_pts);
      }
    }

    // if we dont have any frames with more than min_common_pts we add just the best
    if(vec_frames_pts.empty())
    {
      vec_frames_pts.emplace_back(std::make_pair(max_common_pts, bestFrame));
      // Add connection to new frame in the destination frame
      bestFrame->addFrameVisibilityConnection(this, max_common_pts);

    }

    // Sort all the frames
    const size_t n_frames = vec_frames_pts.size();
    std::sort(vec_frames_pts.begin(),vec_frames_pts.end());

    ordered_connected_frames.resize(n_frames);
    ordered_connected_frames_weights.resize(n_frames);

    for (size_t i = 0; i < n_frames; i++)
    {
      ordered_connected_frames[i] = vec_frames_pts[n_frames-i-1].second;
      ordered_connected_frames_weights[i] = vec_frames_pts[n_frames-i-1].first;
      std::cout<<"FrameA "<<frameId_<<" :: "<<ordered_connected_frames[i]->frameId_<<" :: "<<ordered_connected_frames_weights[i]<<"\n";
    }
  }


  void getFrameVisibilityConnections(std::vector<Frame *> & temp_ordered_connected_frames)
  {
    std::cout<<"GET Frame covisibility: "<<frameId_<<"\n";
    // Loop through all matched points of the frame and vote for each frame that sees any of the points
    Hash_Map<Frame *, size_t> temp_connected_frames_weight;
    for (MapLandmark * map_point: map_points_)
    {
      if (!map_point)
        continue;

      LandmarkObservations & map_obs = map_point->obs_;
      for (auto obs: map_obs)
      {
        if (obs.first == frameId_)
          continue;
        temp_connected_frames_weight[obs.second.frame_ptr]++;
      }
    }
    if (temp_connected_frames_weight.empty())
      return;

    // Add all connections to frames that have at least min_points in common
    int min_common_pts = 15;
    Frame * bestFrame = nullptr;
    size_t max_common_pts = 0;

    std::vector<std::pair<size_t,Frame*> > vec_frames_pts;
    vec_frames_pts.reserve(temp_connected_frames_weight.size());

    for (auto m_iter : temp_connected_frames_weight)
    {
      size_t & n_common_pts = m_iter.second;
      if (n_common_pts > max_common_pts)
      {
        max_common_pts = n_common_pts;
        bestFrame = m_iter.first;
      }

      if (n_common_pts > min_common_pts)
      {
        vec_frames_pts.emplace_back(std::make_pair(n_common_pts, m_iter.first));
        // Add connection to new frame in the destination frame
        m_iter.first->addFrameVisibilityConnection(this, n_common_pts);
      }
    }

    // if we dont have any frames with more than min_common_pts we add just the best
    if(vec_frames_pts.empty())
    {
      vec_frames_pts.emplace_back(std::make_pair(max_common_pts, bestFrame));
      // Add connection to new frame in the destination frame
      //bestFrame->addFrameVisibilityConnection(this, max_common_pts);

    }

    // Sort all the frames
    const size_t n_frames = vec_frames_pts.size();
    std::sort(vec_frames_pts.begin(),vec_frames_pts.end());

    temp_ordered_connected_frames.resize(n_frames);
    //std::vector<size_t> temp_ordered_connected_frames_weights(n_frames);

    for (size_t i = 0; i < n_frames; i++)
    {
      temp_ordered_connected_frames[i] = vec_frames_pts[n_frames-i-1].second;
      //temp_ordered_connected_frames_weights[i] = vec_frames_pts[n_frames-i-1].first;
      //std::cout<<"FrameA "<<frameId_<<" :: "<<temp_ordered_connected_frames[i]->frameId_<<" :: "<<temp_ordered_connected_frames_weights[i]<<"\n";
    }
  }



  void addFrameVisibilityConnection(Frame * frame, size_t & weight)
  {
    connected_frames_weight[frame] = weight;
    updateBestFrameVisibilityConnections();
  }

  void updateBestFrameVisibilityConnections()
  {
    std::cout<<"Update Best covisibility: "<<frameId_<<"\n";
    // Add all connections to frames that have at least min_points in common
    size_t min_common_pts = 15;
    Frame * bestFrame = nullptr;
    size_t max_common_pts = 0;

    std::vector<std::pair<size_t,Frame*> > vec_frames_pts;
    vec_frames_pts.reserve(connected_frames_weight.size());

    for (auto m_iter : connected_frames_weight)
    {
      vec_frames_pts.emplace_back(std::make_pair(m_iter.second, m_iter.first));
    }

    // Sort all the frames
    const size_t n_frames = vec_frames_pts.size();
    std::sort(vec_frames_pts.begin(),vec_frames_pts.end());

    ordered_connected_frames.resize(n_frames);
    ordered_connected_frames_weights.resize(n_frames);

    for (size_t i = 0; i < n_frames; i++)
    {
      ordered_connected_frames[i] = vec_frames_pts[n_frames-i-1].second;
      ordered_connected_frames_weights[i] = vec_frames_pts[n_frames-i-1].first;
      std::cout<<"FrameB "<<frameId_<<" :: "<<ordered_connected_frames[i]->frameId_<<" :: "<<ordered_connected_frames_weights[i]<<"\n";
    }
  }

  std::vector<Frame *> getBestCovisibilityFrames(const size_t &i)
  {
    if (ordered_connected_frames.size() < i)
      return ordered_connected_frames;
    else
      return std::vector<Frame *>(ordered_connected_frames.begin(), ordered_connected_frames.begin()+i);
  }
};


} // namespace VSSLAM
} // namespace openMVG


#endif // FRAME_VSSLAM_HPP}

