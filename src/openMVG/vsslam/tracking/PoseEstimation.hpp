
#pragma once

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/solver_homography_kernel.hpp"

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include <openMVG/numeric/numeric.h>

namespace openMVG  {
namespace VSSLAM  {

using namespace openMVG::robust;
using namespace openMVG::geometry;

static const size_t ACRANSAC_ITER = 64;

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
    ACRansacInfo = ACRANSAC(kernel, vec_inliers, ACRANSAC_ITER, &pE, max_threshold*max_threshold, false);

    if (vec_inliers.size() > SolverType::MINIMUM_SAMPLES *2.5)
    {
      return true;
    }
    return false;
  }

  bool estimate_Rt_fromE
  (
    const Mat3 & K1,
    const Mat3 & K2,
    const Mat & x1,
    const Mat & x2,
    const Mat3 & E,
    std::vector<size_t> & vec_inliers,
    const size_t min_points,
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

  bool estimateRobustRelativePoseHE
  (
    const std::shared_ptr<Frame> & frame_1,
    const std::shared_ptr<Frame> & frame_2,
    Hash_Map<size_t,size_t>  & putative_matches_2_1_idx,
    Mat3 & R,
    Vec3 & t,
    std::vector<size_t> & vec_inliers,
    double & AC_threshold
  )
  {
    // Max error in px for H/E to be considered an inlier
    float max_thresh_px_model = 4;
    // Min points needed to be successfuly triangulated
    size_t min_init_triangulated_pts = 50;


    std::cout<<"Try computing H and E\n";
    double startTimeA = omp_get_wtime();

    // Try to estimate the H and F/E from mathces
    const size_t n_putative_matches = putative_matches_2_1_idx.size();
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
      Hash_Map<size_t,size_t>::iterator m_iter = putative_matches_2_1_idx.begin();
      std::advance(m_iter,m_i);
      pt2D_frame1.col(m_i) = frame_1->pts_undist_[m_iter->second].cast<double>();
      pt2D_frame2.col(m_i) = frame_2->pts_undist_[m_iter->first].cast<double>();
    }

    double stopTimeA = omp_get_wtime();
    double secsElapsedA = stopTimeA - startTimeA; // that's all !
    std::cout<<"Copy data to matrix: "<<secsElapsedA<<"\n";

    // ----------------------------
    // Get camera info
    // ----------------------------
    const Pinhole_Intrinsic * cam_1 = dynamic_cast<const Pinhole_Intrinsic*>(frame_1->cam_intrinsic_);
    const Pinhole_Intrinsic * cam_2 = dynamic_cast<const Pinhole_Intrinsic*>(frame_2->cam_intrinsic_);
    const Mat & cam_1_K = cam_1->K();
    const std::pair<size_t,size_t> cam_1_img(cam_1->w(), cam_1->h());
    const std::pair<size_t,size_t> cam_2_img(cam_2->w(), cam_2->h());
    const Mat & cam_2_K = cam_2->K();


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

    stopTimeA = omp_get_wtime();
    secsElapsedA = stopTimeA - startTimeB; // that's all !
    std::cout<<"Compute H: "<<secsElapsedA<<"\n";

    std::cout<<"Compute H: "<<bValid_H<<" thresh:"<<infoACR_H.first<<" NFA: "<<infoACR_H.second<<" inliers: "<<vec_inliers_H.size()<<"\n";

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

    stopTimeA = omp_get_wtime();
    secsElapsedA = stopTimeA - startTimeB; // that's all !
    std::cout<<"Compute E: "<<secsElapsedA<<"\n";

    std::cout<<"Compute E: "<<bValid_E<<" thresh:"<<infoACR_E.first<<" NFA: "<<infoACR_E.second<<" inliers: "<<vec_inliers_E.size()<<"\n";

    // Check if either of models is ok
    if (!bValid_H && !bValid_E)
    {
      // Neither models were successful! - skip this frame
      std::cout<<"No models available - ABORT initialization and try with next frame\n";
      return false;
    }


    // ----------------------------
    // Estimate Rts
    // ----------------------------
    // Decide which model to use based on NFA (smaller is better)
    bool bUseH = (infoACR_H.second < infoACR_E.second);

    if (bUseH)
    {
      std::cout<<"Use H!\n";
      return false;
    }
    else
    {
      std::cout<<"Use E!\n";
      // Try to get the pose from Essential matrix
      if (estimate_Rt_fromE(cam_1_K, cam_1_K, pt2D_frame1, pt2D_frame2, E, vec_inliers_E, min_init_triangulated_pts, R, t))
      {
        vec_inliers = vec_inliers_E;
        AC_threshold = infoACR_E.first;
        // Successful
        std::cout<<"Motion estimated: "<< R <<"\n"<<"t: "<<t<<"\nc: "<< -R.transpose() * (t)<<"\n";
        return true;
      }
      else
      {
        return false;
      }
    }
  }

  // Find transformation between cameras cam_k and cam_ref: (cam_k)_T_(cam_ref)
  // Look from both directions and find the first intersection
  bool findPoseTransformation
  (
    Frame * cam_k,
    Frame * cam_ref,
    Mat4 & T
  )
  {
    // camera is the reference frame (e.g. transformation is identity)
    if (cam_k == cam_ref)
    {
      T = Mat4::Identity();
      return true;
    }

    // Camera K is already expressed in cam reference
    if (cam_k && cam_k->ref_frame_ == cam_ref)
    {
      T = cam_k->pose_.transformation();
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
            break;

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
        T_backwards_list.push_back(cam_bckw);
        cam_bckw = cam_bckw->ref_frame_;
      }

      // Both cameras come to world origin we just need to traverse both
      if (!cam_frwd && !cam_bckw)
      {
        frwd_match_iter = T_forwards_list.end();
        bcwd_match_iter = T_backwards_list.rend();
      }
    }

    T = Mat4::Identity();
    // Traverse forward list up to where we found camera that is alre
    std::deque<Frame *>::iterator cam_f_iter = T_forwards_list.begin();
    while (cam_f_iter != frwd_match_iter)
    {
      T = T * (*cam_f_iter)->pose_.transformation();
      cam_f_iter++;
    }

    while ((bcwd_match_iter++) != T_backwards_list.rend())
    {
      T = T * (*bcwd_match_iter)->pose_.transformation().inverse();
    }

    return true;
  }

  // Expresses point represented in ref_point_frame in new_point_frame
  bool getRelativePointPosition(Vec3 &pt, Frame * ref_point_frame, Vec3 &pt_new, Frame * new_point_frame )
  {
    if (ref_point_frame == new_point_frame)
    {
      pt_new = pt;
      return true;
    }
    Mat4 T;
    if (findPoseTransformation(new_point_frame,ref_point_frame,T))
    {
      pt_new = Vec4(T * pt.homogeneous()).hnormalized();
      return true;
    }

    // Something went wrong
    return false;

  }
}
}
