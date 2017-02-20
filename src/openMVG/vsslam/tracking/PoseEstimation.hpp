
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/solver_homography_kernel.hpp"
#include <openMVG/vsslam/Camera.hpp>

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/solver_resection_kernel.hpp"

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include <openMVG/numeric/numeric.h>
#include <openMVG/vsslam/Camera.hpp>


namespace openMVG  {
namespace VSSLAM  {

using namespace openMVG::robust;
using namespace openMVG::geometry;

static const size_t ACRANSAC_ITER = 64;

struct ResectionSquaredResidualError {
  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  // Return the squared error
  static double Error(const Mat34 & P, const Vec2 & pt2D, const Vec3 & pt3D){
    return (Project(P, pt3D) - pt2D).squaredNorm();
  }
};

  static double computeHomographyScore
  (
      const Mat & model,
      const Mat & x1,
      const Mat & x2,
      const double gamma
  )
  {
    Mat model_inv = model.inverse();

    double score = 0.0;
    /*
    Mat d2_1in2 = (x2 - ((model * x1.colwise().homogeneous()).colwise().hnormalized())).colwise().squaredNorm();
    Mat d2_2in1 = (x1 - ((model_inv * x2.colwise().homogeneous()).colwise().hnormalized())).colwise().squaredNorm();

    for (size_t pt_i = 0; pt_i < x1.cols(); ++pt_i)
    {
      const double d2_1_2 = d2_1in2(pt_i);
      const double d2_2_1 = d2_2in1(pt_i);
      if (d2_1_2 < gamma)
        score += d2_1_2;
      if (d2_2_1 < gamma)
        score += d2_2_1;
    }*/

    for (size_t pt_i = 0; pt_i < x1.cols(); ++pt_i)
    {
      const Vec2 & pt_1 = x1.col(pt_i);
      const Vec2 & pt_2 = x2.col(pt_i);
      // x2 = H * x1
      Vec2 proj_x = Vec3(model * pt_1.homogeneous()).hnormalized();
      double d2 = (pt_2 - proj_x).squaredNorm();

      if (d2 < gamma)
        score += (gamma - d2);

      // x1 = H^-1 * x2 hnormalized()
      proj_x = Vec3(model_inv * pt_2.homogeneous()).hnormalized();
      d2 = (pt_1 - proj_x).squaredNorm();

      if (d2 < gamma)
        score += (gamma - d2);

    }

    return score;
  }

  static double computeEpipolarScore
(
    const Mat & model,
    const Mat & x1,
    const Mat & x2,
    const double gamma
)
{
  Mat model_inv = model.inverse();

  double score = 0.0;

  /*
  Mat FFX_1 = model * x1.colwise().homogeneous();
  Mat d2_1in2 = ((FFX_1.cwiseProduct(x2.colwise().homogeneous())).colwise().sum().array().pow(2.0f)) / (FFX_1.topRows(2).colwise().squaredNorm()).array();
  Mat FFX_2 = model_inv * x2.colwise().homogeneous();
  Mat d2_2in1 = ((FFX_2.cwiseProduct(x1.colwise().homogeneous())).colwise().sum().array().pow(2.0f)) / (FFX_2.topRows(2).colwise().squaredNorm()).array();

  for (size_t pt_i = 0; pt_i < x1.cols(); ++pt_i)
  {
    const double d2_1_2 = d2_1in2(pt_i);
    const double d2_2_1 = d2_2in1(pt_i);
    if (d2_1_2 < gamma)
      score += d2_1_2;
    if (d2_2_1 < gamma)
      score += d2_2_1;
  }*/

  for (size_t pt_i = 0; pt_i < x1.cols(); ++pt_i)
  {
    const Vec3 pt_2 = x2.col(pt_i).homogeneous();
    const Vec3 pt_1 = x1.col(pt_i).homogeneous();

    // x2T * F * x1 = 0
    Vec3 F_x = model * pt_1;
    double dF_x = F_x.dot(pt_2);
    double d2 = (dF_x * dF_x) /  F_x.head<2>().squaredNorm();

    if (d2 < gamma)
      score += (gamma - d2);

    // x1T * F^-1 * x2 = 0
    F_x = model_inv * pt_2;
    dF_x = F_x.dot(pt_1);
    d2 = (dF_x * dF_x) /  F_x.head<2>().squaredNorm();

    if (d2 < gamma)
      score += (gamma - d2);

  }

  return score;
}


  static bool computeH
  (
    const Mat & x1,
    const Mat & x2,
    const std::pair<size_t,size_t> & img1_size,
    const std::pair<size_t,size_t> & img2_size,
    const double max_threshold,
    Mat3 & H_21,  //x2 = H21 * x1
    std::pair<double, double> & ACRansacInfo,
    std::vector<size_t> & vec_inliers
  )
  {
    //-- Homography robust estimation
    using KernelType =
      ACKernelAdaptor<
        openMVG::homography::kernel::FourPointSolver,
        openMVG::homography::kernel::AsymmetricError,
        UnnormalizerI,
        Mat3>;

    KernelType kernel(
      x1, img1_size.first, img1_size.second,
      x2, img2_size.first, img2_size.second,
      false); // configure as point to point error model.

    // ACRansacOut.first (threshold) ACRansacOut.second(NFA score)
    // Precision is squared distance to projection
    ACRansacInfo = ACRANSAC(kernel, vec_inliers, ACRANSAC_ITER, &H_21, max_threshold*max_threshold, false);

    if (vec_inliers.size() > KernelType::MINIMUM_SAMPLES *2.5)
    {
      return true;
    }
    return false;
  }

  static bool computeE
  (
    const Mat3 & K1,
    const Mat3 & K2,
    const Mat & x1,
    const Mat & x2,
    const std::pair<size_t,size_t> & img1_size,
    const std::pair<size_t,size_t> & img2_size,
    const double max_threshold,
    Mat3 & pE,
    std::pair<double, double> & ACRansacInfo,
    std::vector<size_t> & vec_inliers
  )
  {
    // Use the 5 point solver to estimate E
    using SolverType = openMVG::essential::kernel::FivePointKernel;
    // Define the AContrario adaptor
    using KernelType =
      ACKernelAdaptorEssential<
        SolverType,
        openMVG::fundamental::kernel::EpipolarDistanceError,
        Mat3>;

    KernelType kernel(
      x1, img1_size.first, img1_size.second,
      x2, img2_size.first, img2_size.second,
      K1, K2);

    // Robustly estimation of the Essential matrix and it's precision
    // ACRansacOut.first (threshold) ACRansacOut.second(NFA score)
    // Precision is squared distance to epipolar line
    ACRansacInfo = ACRANSAC(kernel, vec_inliers, ACRANSAC_ITER, &pE, max_threshold*max_threshold, false);

    if (vec_inliers.size() > SolverType::MINIMUM_SAMPLES *2.5)
    {
      return true;
    }
    return false;
  }

  static bool estimate_Rt_from_E
  (
    const Mat3 & K1,
    const Mat3 & K2,
    const Mat & x1,
    const Mat & x2,
    const Mat3 & E,
    std::vector<size_t> & vec_inliers,
    const size_t min_points,
    const double min_cos_angle_points,
    Mat3 & R,
    Vec3 & t
  )
  {
    // Accumulator to find the best solution
    std::vector<size_t> f(4, 0);

    std::vector<Mat3> Rs;  // Rotation matrix.
    std::vector<Vec3> ts;  // Translation matrix.
    std::vector<std::vector<size_t> > pt3Ds_indices(4); // Possible 3D reconstructions

    // Recover best rotation and translation from E.
    openMVG::MotionFromEssential(E, &Rs, &ts);

    //-> Test the 4 solutions will all the points
    assert(Rs.size() == 4);
    assert(ts.size() == 4);

    Mat34 P1;
    Mat3 R1 = Mat3::Identity();
    Vec3 t1 = Vec3::Zero();
    //Vec3 O1 = Vec3::Zero(); // Origin of Cam 1

    P_From_KRt(K1, R1, t1, &P1);

    for (unsigned int i = 0; i < 4; ++i)
    {
      pt3Ds_indices[i].reserve(vec_inliers.size());
      const Mat3 &R2 = Rs[i];
      const Vec3 &t2 = ts[i];
      Mat34 P2;
      P_From_KRt(K2, R2, t2, &P2);

      // 3D point
      Vec3 X;

      for (size_t k = 0; k < vec_inliers.size(); ++k)
      {
        const Vec2
          & x1_ = x1.col(vec_inliers[k]),
          & x2_ = x2.col(vec_inliers[k]);
        TriangulateDLT(P1, x1_, P2, x2_, &X);
        // Test if point is front to the two cameras.
        if (Depth(R1, t1, X) < 0 || Depth(R2, t2, X) < 0)
        {
          continue;
        }

        // Compute ray angle between points
        const Vec3 ray1 = Vec3(R1.transpose() * Vec3( K1.inverse() * Vec3( x1_( 0 ), x1_( 1 ), 1.0 ) )).normalized();
        const Vec3 ray2 = Vec3(R2.transpose() * Vec3( K2.inverse() * Vec3( x2_( 0 ), x2_( 1 ), 1.0 ) )).normalized();
        const double mag = ray1.norm() * ray2.norm();
        const double dotAngle = ray1.dot( ray2 );
        const double cosParallax = clamp( dotAngle / mag, -1.0 + 1.e-8, 1.0 - 1.e-8 );
        // if too small we dont consider it
        if (cosParallax > min_cos_angle_points)
          continue;

        pt3Ds_indices[i].emplace_back(vec_inliers[k]);
        ++f[i];
      }
    }
    // Check the solution:
    const std::vector<size_t>::iterator iter_max = std::max_element(f.begin(), f.end());
    if (*iter_max == 0 || *iter_max < min_points)
    {
      // There is no right solution with points in front of the cameras
      return false;
    }

    // Check that other solutions are not close to the best -> best has to be clear
    const size_t limit_similar = *iter_max * 0.7;

    for (std::vector<size_t>::iterator iter = f.begin(); iter != f.end(); ++iter)
    {
      if (iter != iter_max && (*iter) >limit_similar)
        return false;
    }

    // Return best solution
    const size_t index = std::distance(f.begin(), iter_max);
    R = Rs[index];
    t = ts[index];
    vec_inliers = pt3Ds_indices[index];

    return true;
  }

  static bool estimateRobustRelativePoseHE
  (
    Frame * frame_1,
    Frame * frame_2,
    Hash_Map<size_t,size_t>  & putative_matches_1_2_idx,
    const double min_cos_angle_init_points,
    const size_t min_init_triangulated_pts,
    Mat3 & R,
    Vec3 & t,
    std::vector<size_t> & vec_inliers,
    double & AC_threshold
  )
  {
    // Max error in px for H/E to be considered an inlier
    float max_thresh_px_model = 4;

    std::cout<<"Pose: Try estimating H and E\n";

    double startTimeA = omp_get_wtime();
    // ----------------------------
    // Get camera info
    // ----------------------------
    Camera * cam_1 = frame_1->cam_;
    Camera * cam_2 = frame_2->cam_;
    const Pinhole_Intrinsic * cam_intrinsic_1 = dynamic_cast<const Pinhole_Intrinsic*>(cam_1->cam_intrinsic_ptr);
    const Pinhole_Intrinsic * cam_intrinsic_2 = dynamic_cast<const Pinhole_Intrinsic*>(cam_2->cam_intrinsic_ptr);
    const Mat & cam_1_K = cam_intrinsic_1->K();
    const std::pair<size_t,size_t> cam_1_img(cam_intrinsic_1->w(), cam_intrinsic_1->h());
    const std::pair<size_t,size_t> cam_2_img(cam_intrinsic_2->w(), cam_intrinsic_2->h());
    const Mat & cam_2_K = cam_intrinsic_2->K();

    // Try to estimate the H and F/E from mathces
    const size_t n_putative_matches = putative_matches_1_2_idx.size();
    Mat2X pt2D_frame1(2, n_putative_matches);
    Mat2X pt2D_frame2(2, n_putative_matches);

    // ----------------------------
    // Copy data of matches
    // ----------------------------
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int m_i = 0; m_i < n_putative_matches; ++m_i)
    {
      Hash_Map<size_t,size_t>::iterator m_iter = putative_matches_1_2_idx.begin();
      std::advance(m_iter,m_i);
      pt2D_frame1.col(m_i) = (cam_1->bCalibrated ? frame_1->getFeaturePosition(m_iter->first) : cam_intrinsic_1->remove_disto(frame_1->getFeaturePosition(m_iter->first))).cast<double>();
      pt2D_frame2.col(m_i) = (cam_2->bCalibrated ? frame_2->getFeaturePosition(m_iter->second) : cam_intrinsic_2->remove_disto(frame_2->getFeaturePosition(m_iter->second))).cast<double>();
    }

    // ----------------------------
    // Estimate H
    // ----------------------------
    double startTimeB = omp_get_wtime();

    Mat3 H;
    std::vector<size_t> vec_inliers_H;
    std::pair<double, double> infoACR_H; // first - threshold;  second - NFA score

    bool bValid_H = computeH(pt2D_frame1,pt2D_frame2, cam_1_img, cam_2_img, max_thresh_px_model, H, infoACR_H, vec_inliers_H);
    // Check if enough inliers
    if (bValid_H && vec_inliers_H.size() < min_init_triangulated_pts)
      bValid_H = false;

    double stopTimeB = omp_get_wtime();
    double secsElapsedB = stopTimeB - startTimeB; // that's all !
    std::cout<<"Pose: Computing H ("<<secsElapsedB<<")\n";

    std::cout<<"Pose: H: "<<bValid_H<<" thresh:"<<infoACR_H.first<<" NFA: "<<infoACR_H.second<<" inliers: "<<vec_inliers_H.size()<<"\n";

    // ----------------------------
    // Estimate  E
    // ----------------------------
    startTimeB = omp_get_wtime();

    Mat3 E;
    std::vector<size_t> vec_inliers_E;
    std::pair<double, double> infoACR_E; // first - threshold;  second - NFA score

    bool bValid_E = computeE(cam_1_K, cam_2_K, pt2D_frame1, pt2D_frame2, cam_1_img, cam_2_img, max_thresh_px_model, E, infoACR_E, vec_inliers_E);
    // Check if enough inliers
    if (bValid_E && vec_inliers_E.size() < min_init_triangulated_pts)
      bValid_E = false;

    stopTimeB = omp_get_wtime();
    secsElapsedB = stopTimeB - startTimeB; // that's all !
    std::cout<<"Pose: Computing E ("<<secsElapsedB<<")\n";

    std::cout<<"Pose: E: "<<bValid_E<<" thresh:"<<infoACR_E.first<<" NFA: "<<infoACR_E.second<<" inliers: "<<vec_inliers_E.size()<<"\n";

    // Check if either of models is ok
    if (!bValid_H && !bValid_E)
    {
      // Neither models were successful! - skip this frame
      std::cout<<"Pose: No models available - ABORT initialization and try with next frame\n";
      return false;
    }


    // ----------------------------
    // Estimate Rts
    // ----------------------------
    // Decide which model to use based on NFA (smaller is better)
    bool bUseH = (infoACR_H.second < infoACR_E.second);

    if (bUseH)
    {
      std::cout<<"Pose: Use H!\n";
      return false;
    }
    else
    {
      std::cout<<"Pose: Use E!\n";
      // Try to get the pose from Essential matrix
      if (VSSLAM::estimate_Rt_from_E(cam_1_K, cam_1_K, pt2D_frame1, pt2D_frame2, E, vec_inliers_E, min_init_triangulated_pts, min_cos_angle_init_points, R, t))
      {
        vec_inliers = vec_inliers_E;
        AC_threshold = infoACR_E.first;
        std::cout<<"Pose: Motion estimated OK!\nR:"<<R<<"\nt:"<<t<<"\n";
        return true;
      }
      else
      {
        return false;
      }
    }
  }

  static bool estimateRobustPose
  (
    Frame * frame,
    Hash_Map<MapLandmark* ,size_t> & putative_matches_3D_ptr_cur_idx,
    const size_t min_points,
    Mat3 & R,
    Vec3 & t,
    std::vector<size_t> & vec_inliers,
    double & AC_threshold
  )
  {
    // Max error in px for H/E to be considered an inlier
    float max_thresh_px_model = 4;

    double startTimeA = omp_get_wtime();

    // Try to estimate the H and F/E from mathces
    const size_t n_putative_matches = putative_matches_3D_ptr_cur_idx.size();

    Mat2X pt2D_frame(2, n_putative_matches);
    Mat3X pt3D_frame(3, n_putative_matches);

    const Pinhole_Intrinsic * cam_intrinsic = dynamic_cast<const Pinhole_Intrinsic*>(frame->getCameraIntrinsics());
    const bool bCamCalibrated = frame->getCamCalibrated();
    const Mat3 frame_K_inv = frame->getK_inv();

    // ----------------------------
    // Copy data of matches
    // ----------------------------
    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int m_i = 0; m_i < n_putative_matches; ++m_i)
    {
      Hash_Map<MapLandmark*,size_t>::iterator m_iter = putative_matches_3D_ptr_cur_idx.begin();
      std::advance(m_iter,m_i);
      MapLandmark * map_point = m_iter->first;

      pt2D_frame.col(m_i) = (bCamCalibrated ? frame->getFeaturePosition(m_iter->second) : cam_intrinsic->remove_disto(frame->getFeaturePosition(m_iter->second))).cast<double>();
      pt3D_frame.col(m_i) = (map_point->X_).cast<double>();
    }


    // EPnP

    // Since K calibration matrix is known, compute only [R|t]
    using KernelType =
      openMVG::robust::ACKernelAdaptorResection_K<
      openMVG::euclidean_resection::kernel::EpnpSolver,
        ResectionSquaredResidualError,
        Mat34>;

    KernelType kernel(pt2D_frame, pt3D_frame, cam_intrinsic->K());
    Mat34 P;
    const std::pair<double,double> ACRansacOut =
    openMVG::robust::ACRANSAC(kernel, vec_inliers,1024, &P, max_thresh_px_model*max_thresh_px_model, true);
     // Update the upper bound precision of the model found by AC-RANSAC
     AC_threshold = ACRansacOut.first;

     std::cout<<"PnP inliers: "<<vec_inliers.size()<<"\n";
     const bool bResection = (vec_inliers.size() > 2.5 * openMVG::euclidean_resection::kernel::EpnpSolver::MINIMUM_SAMPLES);

     if (bResection)
     {
       Mat3 K_t;
       KRt_From_P(P, &K_t, &R, &t);
       std::cout<<"Resection: K:"<<K_t<<" \n R:"<<R<<"\n t:"<<t<<"\n";
       return true;
     }
     std::cout<<"Resection: failed!: "<<vec_inliers.size()<<"\n";
     return false;

    /*//P3P
    //--
    // Since K calibration matrix is known, compute only [R|t]
    using KernelType =
      openMVG::robust::ACKernelAdaptorResection_K<
      openMVG::euclidean_resection::P3PSolver,
        ResectionSquaredResidualError,
        Mat34>;

    KernelType kernel(pt2D_frame, pt3D_frame, cam_intrinsic->K());
    Mat34 P;

    // Robust estimation of the Projection matrix and it's precision
    // Precision is squared distance to epipolar line
    const std::pair<double,double> ACRansacOut =
      openMVG::robust::ACRANSAC(kernel, vec_inliers,1024, &P, max_thresh_px_model*max_thresh_px_model, true);
    // Update the upper bound precision of the model found by AC-RANSAC
    AC_threshold = ACRansacOut.first*ACRansacOut.first;

    std::cout<<"P3P inliers: "<<vec_inliers.size()<<"\n";

    // Test if the mode support some points (more than those required for estimation)
    const bool bResection = (vec_inliers.size() > 2.5 * openMVG::euclidean_resection::P3PSolver::MINIMUM_SAMPLES);

    if (bResection)
    {
      Mat3 K_t;
      KRt_From_P(P, &K_t, &R, &t);
      std::cout<<"Resection: K:"<<K_t<<" \n R:"<<R<<"\n t:"<<t<<"\n";
      return true;
    }
    std::cout<<"Resection: failed!: "<<vec_inliers.size()<<"\n";
    return false;
    */

}

  // Find transformation between cameras cam_k and cam_ref: (cam_k)_T_(cam_ref)

  // Look from both directions and find the first intersection
  static bool getRelativeCameraTransformation
  (
    Frame * cam_k,
    Frame * cam_ref,
    Mat4 & T
  )
  {
    // camera is the reference frame (e.g. transformation is identity)
    if (cam_k == cam_ref)
    {
      //std::cout<<"Return identity\n";
      T = Mat4::Identity();
      return true;
    }

    // Camera K is already expressed in cam reference
    if (cam_k && cam_k->ref_frame_ == cam_ref)
    {
      //std::cout<<"Already reference cam\n";
      T = cam_k->getTransformationMatrix_cr();
      return true;
    }

    Frame * cam_frwd = cam_k;  // Find connection between cam_k and cam_ref through T
    Frame * cam_bckw = cam_ref; // Find connection between cam_ref and cam_K through T


    std::deque<Frame *> T_forwards_list;
    std::deque<Frame *>::iterator frwd_match_iter;
    std::deque<Frame *> T_backwards_list;
    std::deque<Frame *>::reverse_iterator bcwd_match_iter;

    while (cam_frwd || cam_bckw)
    {
      if (cam_frwd)
      {
        // Search if we already passed it with reverse path
        bcwd_match_iter = T_backwards_list.rbegin();
        while (bcwd_match_iter != T_backwards_list.rend())
        {
          // We have the reverse path from cam_fwrd in the other list
          // starting with bcwd_match_iter up to beginning of list
          if (cam_frwd == *(bcwd_match_iter))
          {
            bcwd_match_iter++;
            break;
          }
          bcwd_match_iter++;
        }

        // We found a camera in backwards list
        // complete path is:
        //    - complete forward lists (cam_k -> cam_frwd)
        //    - backwards list from (bcwd_match_iter = cam_frwd -> T_backwards_list.begin()=cam_ref) needs to be inverted!!
        if (bcwd_match_iter != T_backwards_list.rend())
        {
          frwd_match_iter = T_forwards_list.end();
          break;
        }
        // Add forward cam to the list
        //std::cout<<"Put to TA: "<<(cam_frwd)->getFrameId()<<"\n";
        T_forwards_list.push_back(cam_frwd);
        cam_frwd = cam_frwd->ref_frame_;
      }

      if (cam_bckw)
      {

        // Search if we already passed it with forward path
        frwd_match_iter = T_forwards_list.begin();
        while (frwd_match_iter != T_forwards_list.end())
        {
          // If we find camera in forward list we have the path from
          // cam_k to cam_bckw in forward list and the rest in backwards
          if (cam_bckw == *(frwd_match_iter))
            break;

          frwd_match_iter++;
        }

        // We found a camera in forward list
        // complete path is:
        //    - complete forward lists (cam_k -> frwd_match_iter = cam_bckw)
        //    - backwards list from (T_backwards_list.end() = cam_bckw -> T_backwards_list.begin()=cam_ref) needs to be inverted!!
        if (frwd_match_iter != T_forwards_list.end())
        {
          bcwd_match_iter = T_backwards_list.rend();
          break;
        }

        // Add backwards cam to the list
        //std::cout<<"Put to TB: "<<(cam_bckw)->getFrameId()<<"\n";
        T_backwards_list.push_back(cam_bckw);
        cam_bckw = cam_bckw->ref_frame_;
      }

      // Both cameras come to world origin we just need to traverse both
      if (!cam_frwd && !cam_bckw)
      {
        //std::cout<<"Both came to W\n";
        frwd_match_iter = T_forwards_list.end();
        bcwd_match_iter = T_backwards_list.rbegin();
      }
      //std::cout<<"End iter\n";
    }

    T = Mat4::Identity();
    // Traverse forward list up to where we found camera that is alre
    std::deque<Frame *>::iterator cam_f_iter = T_forwards_list.begin();
    while (cam_f_iter != frwd_match_iter)
    {
      //std::cout<<"Added TA: "<<(*cam_f_iter)->getFrameId()<<"\n";
      T = T * (*cam_f_iter)->getTransformationMatrix_cr();
      cam_f_iter++;
    }

    while ((bcwd_match_iter) != T_backwards_list.rend())
    {
      //std::cout<<"Added TB: "<<(*bcwd_match_iter)->getFrameId()<<"\n";
      T = T * (*bcwd_match_iter)->getTransformationMatrix_rc();
      bcwd_match_iter++;
    }

    return true;
  }


  // Expess point pt with ref_point_frame in new_point_frame (result in new_pt)
  // For now only with global eculidean representation of points and global represenation of poitns
  static bool getRelativePointPosition(const Vec3 & pt, const Frame * ref_point_frame, Vec3 & pt_new, const Frame * new_point_frame )
  {
    if (ref_point_frame == new_point_frame)
    {
      pt_new = pt;
      return true;
    }
    // PT is in {W}
    if (ref_point_frame == nullptr)
    {
      // Cam is in {W}
      if (new_point_frame->ref_frame_ == nullptr)
      {
        pt_new = Vec4(new_point_frame->getTransformationMatrix_cr() * pt.homogeneous()).hnormalized();
        return true;
      }
    }
    /*
    Mat4 T;
    if (getRelativeCameraTransformation(new_point_frame,ref_point_frame,T))
    {
      pt_new = Vec4(T * pt.homogeneous()).hnormalized();
      return true;
    }
  */
    // Something went wrong
    return false;

  }

  static bool computeFundamentalMatrix(Frame * frame_1, Frame * frame_2, Mat3 & F)
  {
    Mat4 T_2_1;
    // Get T between cameras (2_T_1)
    // Transformation from frame_1 to frame_2
    if (!getRelativeCameraTransformation(frame_2,frame_1,T_2_1))
    {
      // Shouldnt happen!
      return false;
    }

    // Get Rt from T
    const double scale = T_2_1.block(0,0,3,1).norm();
    const Mat3 R = T_2_1.block(0,0,3,3)/scale;
    const Vec3 t = T_2_1.block(0,3,3,1)/scale;
    Mat3 tx;
    tx<<0,-t(2),t(1), t(2),0,-t(0), -t(1),t(0),0;

    // Compute essential matrix
    const Mat3 K_inv_1 = frame_1->getK_inv();
    const Mat3 K_inv_2 = frame_2->getK_inv();

    F = K_inv_2.transpose() * tx * R * K_inv_1;
    return true;

  }


}
}
