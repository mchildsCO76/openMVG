
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FRAME_VSSLAM_HPP
#define FRAME_VSSLAM_HPP

#include <iostream>
#include <memory>

#include <openMVG/types.hpp>
#include <openMVG/numeric/numeric.h>
#include <openMVG/features/features.hpp>
#include <openMVG/geometry/Similarity3.hpp>
#include "openMVG/multiview/projection.hpp"
#include <openMVG/cameras/cameras.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include <openMVG/vsslam/Camera.hpp>
//#include <openMVG/vsslam/tracking/PoseEstimation.hpp>

namespace openMVG  {
namespace VSSLAM  {

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;

struct MapLandmark;

/// Frame
class Frame : public std::enable_shared_from_this<Frame>
{
private:
  // Basic stats
  size_t frame_id_;

  // Is in the map -> Active
  bool active_;

public:
  // Camera
  Camera * cam_;

  /// Detected features
  std::unique_ptr<features::Regions> regions_;
  // Covariance of each detection
  std::vector<Mat> pts_cov_;
  // Distorted/undistorted points - if camera is calibrated (for faster reading)
  std::vector<Vec2> pts_undist_;
  // Number of features detected
  size_t n_feats_ = 0;

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

  double AC_reprojection_thresh_ = 0.0f;
  // Owner
  size_t owner_frame_id_ = std::numeric_limits<size_t>::max(); // undefined
  Frame * ref_frame_ = nullptr;


  Frame
  (
    const size_t fId,
    const size_t camId,
    Camera * cam
  ): frame_id_(fId), cam_(cam)
  {
    active_ = false;
    T_cr_ = T_rc_ = Mat4::Identity();

    P_cr_ = dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->K() * T_cr_.block(0,0,3,4);
    P_rc_ = dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->Kinv() * T_rc_.block(0,0,3,4);

    R_rc_ = Mat3::Identity();
    R_cr_ = Mat3::Identity();
  }

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
  const size_t & getFrameId() const
  {
    return frame_id_;
  }
  const size_t & getCamId() const
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

  const size_t & getNumberOfFeatures() const
  {
    return n_feats_;
  }

  Vec2 getFeaturePositionDetected(const size_t & i) const
  {
    return regions_->GetRegionPosition(i);
  }

  const Vec2 & getFeaturePosition(const size_t & i) const
  {
    // If camera is calibrated we precalculate the undistorted points
    return pts_undist_[i];
  }

  void updateFeaturesData()
  {
    // Update number of features detected
    n_feats_ = regions_->RegionCount();
    // Initialize pointers to map points (to nullptr)
    map_points_ = std::vector<MapLandmark *>(n_feats_,static_cast<VSSLAM::MapLandmark *>(nullptr));

    // Undistort points
    // If calibration is known undistort points (save process time later)
    // Otherwise just save distorted points
    pts_undist_.resize(n_feats_);
    undistortPoints();
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

  bool isPointInFrame(const Vec2 & pt) const
  {
    return cam_->isPointInImageBorders(pt);
  }





  void clearMapPoints()
  {
    map_points_ = std::vector<MapLandmark *>(n_feats_,static_cast<VSSLAM::MapLandmark *>(nullptr));
  }

  void clearMapPoint(const size_t & m_i)
  {
    map_points_[m_i] = static_cast<VSSLAM::MapLandmark *>(nullptr);
  }

  size_t getNumberMapPoints() const
  {
    size_t n_matches = 0;
    for (size_t i = 0; i<map_points_.size(); ++i)
    {
      if (map_points_[i])
        n_matches++;
    }
    return n_matches;
  }

  double getSquaredReprojectionError(const Vec3 & pt_frame, const size_t feat_id) const
  {
    const Vec2 & obs_cur = pts_undist_[feat_id];

    const IntrinsicBase * cam_intrinsic = cam_->cam_intrinsic_ptr;
    // Compute residual error in current frame
    // We add distortion to 3D points if we have it
    const Vec2 pt_cur = cam_intrinsic->cam2ima(cam_intrinsic->have_disto()?cam_intrinsic->add_disto(pt_frame.hnormalized()):pt_frame.hnormalized());

    return (obs_cur - pt_cur).squaredNorm();
  }



  // ----------------------
  // -- Set pose: Transformation from r:reference {W} to c:current
  // ----------------------
  /*void setPose_cr(const geometry::Similarity3 & pose_sim)
  {
    // T = [ sR st ; 0 1]
    T_cr_.block(0,0,3,3) = pose_sim.scale_ * pose_sim.pose_.rotation();
    T_cr_.block(0,3,3,1) = pose_sim.scale_ * pose_sim.pose_.translation();
    updateMatrices();
  }
  void setPose_cr(const geometry::Pose3 & pose)
  {
    // T = [ sR st ; 0 1]
    T_cr_.block(0,0,3,3) = pose.rotation();
    T_cr_.block(0,3,3,1) = pose.translation();
    updateMatrices();
  }
*/

  void setPose_cr_T(const Mat4 & this_T_tmp_ref, Frame * tmp_ref = nullptr)
  {
    // we already have the transformation in the reference frame
    if (tmp_ref == ref_frame_)
    {
      T_cr_ = this_T_tmp_ref;
    }
    /*else
    {
      // We have the transformation from tmp_ref to this_frame (this_T_tmp_ref)
      // What we want is the transformation from ref_frame to this_frame (this_T_ref_frame)
      // this_T_ref_frame = this_T_tmp_ref * tmp_ref_T_ref_frame
      Mat4 tmp_ref_T_ref_frame;
      getRelativeCameraTransformation(tmp_ref,ref_frame_,tmp_ref_T_ref_frame);

      T_cr_ = this_T_tmp_ref * tmp_ref_T_ref_frame;
    }*/
    updateMatrices();
  }
  void setPose_cr_Rts(const Mat3 & R, const Vec3 & t, const double & s, Frame * tmp_ref = nullptr)
  {
    // we already have the transformation in the reference frame
    if (tmp_ref == ref_frame_)
    {
      T_cr_.block(0,0,3,3) = s * R;
      T_cr_.block(0,3,3,1) = s * t;
    }
    /*else
    {
      // We have the transformation from tmp_ref to this_frame (this_T_tmp_ref)
      // What we want is the transformation from ref_frame to this_frame (this_T_ref_frame)
      // this_T_ref_frame = this_T_tmp_ref * tmp_ref_T_ref_frame
      Mat4 tmp_ref_T_ref_frame;
      getRelativeCameraTransformation(tmp_ref,ref_frame_,tmp_ref_T_ref_frame);

      T_cr_.block(0,0,3,3) = s * R;
      T_cr_.block(0,3,3,1) = s * t;
      T_cr_ = T_cr_ * tmp_ref_T_ref_frame;

    }*/
    updateMatrices();
  }

  void getPose_cr_Rts(Mat3 & R, Vec3 & t, double & s, Frame * tmp_ref = nullptr) const
  {
    // we already have the pose in the reference frame
    if (tmp_ref == ref_frame_)
    {
      s = T_cr_.block(0,0,3,1).norm();
      R = T_cr_.block(0,0,3,3)/s;
      t = T_cr_.block(0,3,3,1)/s;
    }
    /*else
    {
      // We have the transformation from ref_frame to this_frame (this_T_ref_frame)
      // What we want is the transformation from tmp_ref to this_frame (this_T_tmp_ref)
      // this_T_tmp_ref = this_T_tmp_ref * tmp_ref_T_ref_frame

      // We compute the transformation of the this camera with frame_ref reference frame
      // this_T_frame_ref
      Mat4 this_T_tmp_ref;
      getRelativeCameraTransformation(this,tmp_ref,this_T_tmp_ref);
      s = this_T_tmp_ref.block(0,0,3,1).norm();
      R = this_T_tmp_ref.block(0,0,3,3)/s;
      t = this_T_tmp_ref.block(0,3,3,1)/s;
    }*/
  }


  // ----------------------
  // -- Set reference frame:
  // -- Pose: Transformation from reference (r) to current (c) frame (c_T_r)
  // ----------------------
  void setReferenceFrame (Frame * new_ref_frame)
  {
    ref_frame_ = new_ref_frame;
  }
/*
  void setReferenceFrame_cr (Frame * new_ref_frame, const geometry::Similarity3 &pose_sim)
  {
    ref_frame_ = new_ref_frame;
    //
    setPose_cr(pose_sim);
  }

  void setReferenceFrame_cr (Frame * new_ref_frame, const geometry::Pose3 &pose)
  {
    ref_frame_ = new_ref_frame;
    setPose_cr(pose);
  }*/

  void setReferenceFrame_cr_Rts
  (
    Frame * new_ref_frame,
    const Mat3 & R,
    const Vec3 & t,
    const double s = 1.0
  )
  {
    ref_frame_ = new_ref_frame;
    setPose_cr_Rts(R,t,s, ref_frame_);
  }


  void updateMatrices()
  {
    // Transformation
    T_rc_ = T_cr_.inverse();
    // Rotation
    double s = T_rc_.block(0,0,3,1).norm();
    R_rc_ = T_rc_.block(0,0,3,3)/s;
    s = T_cr_.block(0,0,3,1).norm();
    R_cr_ = T_cr_.block(0,0,3,3)/s;

    // Projection (only useful if its calibrated camera otherwise it changes)
    if (getCamCalibrated())
    {
      P_cr_ = dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->K() * T_cr_.block(0,0,3,4);
      P_rc_ = dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->Kinv() * T_rc_.block(0,0,3,4);
    }
  }

  Mat34 getProjectionMatrix() const
  {
    if (getCamCalibrated())
    {
      return P_cr_;
    }
    else
    {
      return dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->K() * T_cr_.block(0,0,3,4);
    }
  }

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
  void computeFrameVisibilityConnections()
  {
    std::cout<<"Frame: Update frame covisibility: "<<frame_id_<<"\n";
    // Loop through all matched points of the frame and vote for each frame that sees any of the points
    connected_frames_weight.clear();
    for (MapLandmark * map_point: map_points_)
    {
      if (!map_point)
        continue;
      LandmarkObservations & map_obs = map_point->obs_;

      for (auto & obs: map_obs)
      {
        if (obs.first == frame_id_)
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
      //std::cout<<"FrameA "<<frame_id_<<" :: "<<ordered_connected_frames[i]->frame_id_<<" :: "<<ordered_connected_frames_weights[i]<<"\n";
    }
  }

  // Add new Frame connection and update list
  void addFrameVisibilityConnection(Frame * frame, size_t & weight)
  {
    connected_frames_weight[frame] = weight;
    updateBestFrameVisibilityConnections();
  }

  // Update list of best frame visibility connections
  void updateBestFrameVisibilityConnections()
  {
    //std::cout<<"Update Best covisibility: "<<frame_id_<<"\n";
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
      //std::cout<<"FrameB "<<frame_id_<<" :: "<<ordered_connected_frames[i]->frame_id_<<" :: "<<ordered_connected_frames_weights[i]<<"\n";
    }
  }

  // Return n-best frames in covisibility graph
  std::vector<Frame *> getBestCovisibilityFrames(const size_t &n)
  {
    if (ordered_connected_frames.size() < n)
      return ordered_connected_frames;
    else
      return std::vector<Frame *>(ordered_connected_frames.begin(), ordered_connected_frames.begin()+n);
  }

  // Computes covisibility connecctions on the fly (to prevent screwing graph before putting frame in the system)
  void getFrameVisibilityConnections(std::vector<Frame *> & temp_ordered_connected_frames, const size_t N = std::numeric_limits<size_t>::max()) const
  {
    // Min points two frames have to have in common to be considered a pair
    int min_common_pts = 15;

    // Loop through all matched points of the frame and vote for each frame that sees any of the points
    Hash_Map<Frame *, size_t> temp_connected_frames_weight;
    for (MapLandmark * map_point: map_points_)
    {
      if (!map_point)
        continue;

      LandmarkObservations & map_obs = map_point->obs_;
      for (auto obs: map_obs)
      {
        if (obs.first == frame_id_)
          continue;
        temp_connected_frames_weight[obs.second.frame_ptr]++;
      }
    }
    if (temp_connected_frames_weight.empty())
      return;

    // Add all connections to frames that have at least min_points in common
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
        //m_iter.first->addFrameVisibilityConnection(this, n_common_pts);
      }
    }

    // if we dont have any frames with more than min_common_pts we add just the best
    if(vec_frames_pts.empty())
    {
      vec_frames_pts.emplace_back(std::make_pair(max_common_pts, bestFrame));
      // Add connection to new frame in the destination frame
      //bestFrame->addFrameVisibilityConnection(this, max_common_pts);

    }

    // Determine how many frames we want to get back
    const size_t n_frames = vec_frames_pts.size();
    const size_t n_return_frames = (N==0 ? n_frames : std::min<size_t>(n_frames,N));

    // Sort all the frames
    std::sort(vec_frames_pts.begin(),vec_frames_pts.end());

    // Return N best
    temp_ordered_connected_frames.resize(n_return_frames);

    for (size_t i = 0; i < n_return_frames; i++)
    {
      temp_ordered_connected_frames[i] = vec_frames_pts[n_frames-i-1].second;
    }
  }
};


} // namespace VSSLAM
} // namespace openMVG


#endif // FRAME_VSSLAM_HPP}

