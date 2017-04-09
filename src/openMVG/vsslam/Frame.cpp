
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <openMVG/vsslam/Frame.hpp>

namespace openMVG  {
namespace VSSLAM  {

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;


  Frame::Frame
  (
    const IndexT fId,
    const IndexT camId,
    Camera * cam
  ): frame_id_(fId), cam_(cam)
  {
    active_ = false;
    T_cr_ = T_rc_ = Mat4::Identity();

    R_rc_ = R_cr_ = Mat3::Identity();

    P_cr_ = dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->K() * T_cr_.block(0,0,3,4);
    P_rc_ = dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->Kinv() * T_rc_.block(0,0,3,4);


    timestamp_ = fId;

  }

  void Frame::updateFeaturesData()
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

  void Frame::undistortPoints()
  {
    // If we have distortion and can model it we remove it for faster processing later
    if (cam_->cam_intrinsic_->have_disto() && cam_->bCalibrated)
    {
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for
      #endif
      for ( IndexT i = 0; i < n_feats_; ++i )
      {
        pts_undist_[i] = cam_->cam_intrinsic_->remove_disto(regions_->GetRegionPosition(i));
      }
    }
    else
    {
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for
      #endif
      for ( IndexT i = 0; i < n_feats_; ++i )
      {
        pts_undist_[i] = regions_->GetRegionPosition(i);
      }
    }
  }


  bool Frame::isPointInFrame(const Vec2 & pt) const
  {
    return cam_->isPointInImageBorders(pt);
  }

  void Frame::clearMapPoints()
  {
    map_points_ = std::vector<MapLandmark *>(n_feats_,static_cast<VSSLAM::MapLandmark *>(nullptr));
  }

  void Frame::clearMapPoint(const IndexT & p_i)
  {
    map_points_[p_i] = static_cast<VSSLAM::MapLandmark *>(nullptr);
  }

  void Frame::setMapPoint(const IndexT p_i, MapLandmark * ml)
  {
    if (ml && p_i < n_feats_)
    {
      map_points_[p_i] = ml;
      ml->last_local_map_frame_id_ = frame_id_;
    }
  }

  size_t Frame::getNumberMapPoints() const
  {
    IndexT n_matches = 0;
    for (IndexT i = 0; i<map_points_.size(); ++i)
    {
      if (map_points_[i])
        n_matches++;
    }
    return n_matches;
  }


  size_t Frame::getNumberGlobalMapPoints() const
  {
    IndexT n_matches = 0;
    for (IndexT i = 0; i<map_points_.size(); ++i)
    {
      if (map_points_[i])
        n_matches++;
    }
    return n_matches;
  }


  double Frame::getSquaredReprojectionError(const Vec3 & pt_frame, const IndexT feat_id) const
  {
    const Vec2 & obs_cur = pts_undist_[feat_id];

    const IntrinsicBase * cam_intrinsic = cam_->cam_intrinsic_ptr;
    // Compute residual error in current frame
    // We add distortion to 3D points if we have it
    const Vec2 pt_cur = cam_intrinsic->cam2ima(cam_intrinsic->have_disto()?cam_intrinsic->add_disto(pt_frame.hnormalized()):pt_frame.hnormalized());

    return (obs_cur - pt_cur).squaredNorm();
  }

  bool Frame::checkFeatureAssociation(const Vec3 & pt_3D_frame, const IndexT feat_id, double thresh) const
  {
    if (pt_3D_frame(2) < 0)
      return false;

    // Project the point to image plane of frame 2
    const IntrinsicBase * cam_intrinsic = cam_->cam_intrinsic_ptr;
    const Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(cam_intrinsic->have_disto()?cam_intrinsic->add_disto(pt_3D_frame.hnormalized()):pt_3D_frame.hnormalized());

    // Check if projection is actually in the image borders
    if (!isPointInFrame(pt_3D_frame_projected))
      return false;

    // Error between projection and measurement
    const Vec2 pt_err = pts_undist_[feat_id] - pt_3D_frame_projected;

    // Reprojection error
    //return (pt_err.squaredNorm() < thresh);

    //// Chi2 error
    return (pt_err.transpose() * pts_information_mat_[feat_id].cwiseProduct(pts_information_mat_[feat_id]) * pt_err) < thresh;
  }

  void Frame::computeSceneStatistics()
  {
    std::vector<float> vec_depths;
    vec_depths.reserve(n_feats_);

    for (size_t p_i=0; p_i < map_points_.size(); ++p_i)
    {
      if (!map_points_[p_i])
        continue;

      MapLandmark * pt_3D = map_points_[p_i];
      const Vec3 & pt_3D_w = pt_3D->getWorldPosition();
      // Compute z
      float z = (T_cr_.block(2,0,1,4)*Vec4(pt_3D_w.homogeneous())).norm();
      vec_depths.push_back(z);
    }
    std::sort(vec_depths.begin(), vec_depths.end());

    f_scene_median = vec_depths[(vec_depths.size()-1)/2];
    f_scene_095 = vec_depths[((vec_depths.size()-1)*0.85)];

    std::cout<<"F: "<<frame_id_<<"Scene statistics: m: "<<f_scene_median<<" 95: "<<f_scene_095<<" L: "<<vec_depths[vec_depths.size()-1]<<"\n";
  }

  bool Frame::checkLandmarkPosition(const Vec3 & map_landmark_3D_in_frame)
  {
    if (map_landmark_3D_in_frame(2)<0 || map_landmark_3D_in_frame(2) > f_scene_095*3)
      return false;
    return true;

  }


  // Update all pose related matrices
  // Everything is updated from T_cr_
  void Frame::updatePoseMatrices()
  {
    // Transformation [sR t]
    T_rc_ = T_cr_.inverse();

    // Scale
    double s_cr = T_cr_.block(0,0,3,1).norm();
    // Rotation
    R_cr_ = T_cr_.block(0,0,3,3)/s_cr;
    R_rc_ = R_cr_.transpose();

    // Compute camera center in world coordinates
    if (ref_frame_ == nullptr)
    {
      // O_w_ = -1/s_cr * R_cr_.transpose() * t_cr
      Vec3 t_cr = T_cr_.block(0,3,3,1);
      O_w_ = - 1/s_cr * R_rc_ * t_cr;
    }
    else
    {
      // Get poseInverse with world reference frame (nullptr)
      // Translation is actually center of camera in world coordinates
      double s; Mat3 R;
      getPoseInverse_Rts(R,O_w_,s,nullptr);
    }

    // Projection (only useful if its calibrated camera otherwise it changes)
    if (getCamCalibrated())
    {

      double s = T_cr_.block(0,0,3,1).norm();
      Vec3 t = T_cr_.block(0,3,3,1);
      //P_cr_ = dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->K() * T_cr_.block(0,0,3,4);
      P_cr_ = dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->K() * T_cr_.block(0,0,3,4);


      P_rc_ = dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->Kinv() * T_rc_.block(0,0,3,4);
    }
  }

  // Get pose of camera in desired (tmp_ref) reference frame
  // aka get transformation from desired reference (tmp_ref) frame to this frame
  // (e.g. from world to camera reference frame)
  void Frame::getPose_Rts
  (
    Mat3 & R,
    Vec3 & t,
    double & s,
    Frame * tmp_ref
  ) const
  {
    // we already have the pose in the reference frame
    if (tmp_ref == ref_frame_)
    {
      // T = [sR t]
      s = T_cr_.block(0,0,3,1).norm();
      R = R_cr_;
      t = T_cr_.block(0,3,3,1);
    }
    // TODO: relative poses
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
      t = this_T_tmp_ref.block(0,3,3,1);
    }*/
  }

  // Get transformation from this frame to desired reference (tmp_ref) frame
  // (e.g. from camera to world reference frame)
  void Frame::getPoseInverse_Rts
  (
    Mat3 & R,
    Vec3 & t,
    double & s,
    Frame * tmp_ref
  ) const
  {
    // we already have the pose in the reference frame
    if (tmp_ref == ref_frame_)
    {
      // T = [sR t]
      s = T_rc_.block(0,0,3,1).norm();
      R = R_rc_;
      t = T_rc_.block(0,3,3,1);
    }
    // TODO: relative poses
  }

  // Get transformation from this frame to desired reference (tmp_ref) frame expressed in sim3
  // (e.g. from camera to world reference frame)
  void Frame::getPoseInverse_sim3
  (
    Eigen::Matrix<double, 7, 1> & v_state,
    Frame * tmp_ref // Pose [sR t] expressed in tmp_ref reference frame
  ) const
  {
    //double s; Mat3 R; Vec3 t;
    // Get [sR t]
    //getPoseInverse_Rts(R,t,s,tmp_ref);

    // Express in sim3
    Sim3_log(T_rc_,v_state);
  }

  // Get transformation from this frame to desired reference (tmp_ref) frame expressed with state vector
  // (e.g. from camera to world reference frame)
  void Frame::getPoseInverse_StateVector
  (
    Eigen::Matrix<double, 12, 1> & v_state,
    Frame * tmp_ref // Pose [sR t] expressed in tmp_ref reference frame
  ) const
  {
    Eigen::Matrix<double,7,1> log_state;
    // Get sim3 representation of poseInverse
    getPoseInverse_sim3(log_state,tmp_ref);
    // Add it to state vector
    v_state.template head<7>() = log_state;
    // Augment state vector with camera parameters
    const Mat3 K = dynamic_cast<const Pinhole_Intrinsic*>(getCameraIntrinsics())->K();
    // fx,fy,ppx,ppy,d
    // TODO: Check how is with distortion parameters
    v_state.template tail<5>() << K(0,0), K(0,0), K(0,2) , K(1,2), 0.0;
  }

  // Get camera center expressed in world coordinates (regardless of representation of camera)
  const Vec3 & Frame::getCameraCenter() const
  {
    return O_w_;
  }

  // Set pose of the frame with transformation T expressed
  // as transformation from tmp_reference to this frame (camera) : c_T_tmp_ref
  // If reference frame of the given transformation (tmp_ref) is not equal to reference frame of the frame (ref_frame)
  // we need transformation from reference of the frame to reference of given transformation: tmp_ref_T_ref
  // Then we can compute pose: c_T_ref = c_T_tmp_ref * tmp_ref_T_ref
  void Frame::setPose_T
  (
    const Mat4 & c_T_tmp_ref,
    Frame * tmp_ref /* = nullptr */
  )
  {
    // we already have the transformation in the reference frame
    if (tmp_ref == ref_frame_)
    {
      T_cr_ = c_T_tmp_ref;
    }
    /*else
    {
      // If reference frame of the given transformation (tmp_ref) is not equal to reference frame of the frame (ref_frame)
      // we need transformation from reference of the frame to reference of given transformation: tmp_ref_T_ref
      // Then we can compute pose: c_T_ref = c_T_tmp_ref * tmp_ref_T_ref
      Mat4 tmp_ref_T_ref_frame;
      getRelativeCameraTransformation(tmp_ref,ref_frame_,tmp_ref_T_ref_frame);

      T_cr_ = this_T_tmp_ref * tmp_ref_T_ref_frame;
    }*/
    updatePoseMatrices();
  }

  void Frame::setPose_Rts
  (
    const Mat3 & R,
    const Vec3 & t,
    const double & s,
    Frame * tmp_ref /* = nullptr*/
  )
  {
    // we already have the transformation in the reference frame
    if (tmp_ref == ref_frame_)
    {
      T_cr_ = Mat4::Identity();
      T_cr_.block(0,0,3,3) = s * R;
      T_cr_.block(0,3,3,1) = t;
    }
    /*else
    {
      // If reference frame of the given transformation (tmp_ref) is not equal to reference frame of the frame (ref_frame)
      // we need transformation from reference of the frame to reference of given transformation: tmp_ref_T_ref
      // Then we can compute pose: c_T_ref = c_T_tmp_ref * tmp_ref_T_ref

      Mat4 tmp_ref_T_ref_frame;
      getRelativeCameraTransformation(tmp_ref,ref_frame_,tmp_ref_T_ref_frame);

      T_cr_.block(0,0,3,3) = s * R;
      T_cr_.block(0,3,3,1) = t;
      T_cr_ = T_cr_ * tmp_ref_T_ref_frame;

    }*/
    updatePoseMatrices();
  }



  // Set pose of the frame with [sR t] transformation expressed
  // as sim3 vector [upsilon, omega, sigma] from this frame (camera) to reference frame : tmp_ref_T_c
  // TODO: If reference frame of the given transformation (tmp_ref) is not equal to reference frame of the frame (ref_frame)
  void Frame::setPoseInverse_sim3
  (
    const Eigen::Matrix<double, 7, 1> & v_state,
    Frame * tmp_ref
  )
  {
    Mat4 T;
    // From sim3 to Sim3
    Sim3_exp(v_state,T);

    // we already have the transformation in the reference frame
    if (tmp_ref == ref_frame_)
    {
      // Update transformation matrix
      T_cr_ = T.inverse();

    }
    /*
    else
    {

    }
    */
    // Update the rest of the data
    updatePoseMatrices();
  }




  // ----------------------
  // -- Set reference frame:
  // -- Pose: Transformation from reference (r) to current (c) frame (c_T_r)
  // ----------------------
  // Set frame as the reference frame
  void Frame::setReferenceFrame (Frame * new_ref_frame)
  {
    ref_frame_ = new_ref_frame;
  }

  void Frame::setReferenceFrameAndPose_T (Frame * new_ref_frame, const Mat4 &T)
  {
    ref_frame_ = new_ref_frame;
    setPose_T(T,new_ref_frame);
  }

  void Frame::setReferenceFrameAndPose_Rts
  (
    Frame * new_ref_frame,
    const Mat3 & R,
    const Vec3 & t,
    const double s /* = 1.0 */
  )
  {
    ref_frame_ = new_ref_frame;
    setPose_Rts(R,t,s, ref_frame_);
  }




  Mat34 Frame::getProjectionMatrix() const
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


  // ----------------------
  // -- Covisibility stuff
  // ----------------------

  // Compute Frame visibility connections
  void Frame::computeFrameVisibilityConnections()
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

    std::vector<std::pair<IndexT,Frame*> > vec_frames_pts;
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
    }
  }

  // Add new Frame connection and update list
  void Frame::addFrameVisibilityConnection(Frame * frame, size_t & weight)
  {
    connected_frames_weight[frame] = weight;
    updateBestFrameVisibilityConnections();
  }

  // Update list of best frame visibility connections
  void Frame::updateBestFrameVisibilityConnections()
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
  std::vector<Frame *> Frame::getBestCovisibilityFrames(const size_t &n)
  {
    if (ordered_connected_frames.size() < n)
      return ordered_connected_frames;
    else
      return std::vector<Frame *>(ordered_connected_frames.begin(), ordered_connected_frames.begin()+n);
  }

  // Computes covisibility connecctions on the fly (to prevent screwing graph before putting frame in the system)
  void Frame::getFrameVisibilityConnections(std::vector<Frame *> & temp_ordered_connected_frames, const size_t N) const
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

} // namespace VSSLAM
} // namespace openMVG



