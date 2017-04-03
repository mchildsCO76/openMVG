
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
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include <openMVG/numeric/numeric.h>
#include <openMVG/matching/indMatch.hpp>

#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/Camera.hpp>

using namespace openMVG;
using namespace openMVG::robust;
using namespace openMVG::geometry;

namespace openMVG  {
namespace VSSLAM  {

class Camera;
class Frame;

struct ResectionSquaredResidualError {
  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  // Return the squared error
  static double Error(const Mat34 & P, const Vec2 & pt2D, const Vec3 & pt3D){
    return (Project(P, pt3D) - pt2D).squaredNorm();
  }
};

class PoseEstimator
{
public:
  static const size_t ACRANSAC_ITER = 256;

  static double computeHomographyScore
  (
    const Mat & model,
    const Mat & x1,
    const Mat & x2,
    const double gamma
  );

  static double computeEpipolarScore
  (
    const Mat & model,
    const Mat & x1,
    const Mat & x2,
    const double gamma
  );


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
  );

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
  );

  static bool estimate_Rt_from_E
  (
    const Mat3 & K1,
    const Mat3 & K2,
    const Mat & x1,
    const Mat & x2,
    const Mat3 & E,
    std::vector<size_t> & vec_inliers,
    const double min_cos_angle_points,
    Mat3 & R,
    Vec3 & t
  );

  static bool estimateRobustRelativePosePinholeHE
  (
    const Frame * frame_1,
    const Frame * frame_2,
    const openMVG::matching::IndMatches & vec_putative_matches_1_2_idx,
    const double min_cos_angle_init_points,
    const size_t min_init_triangulated_pts,
    Mat4 & T,
    std::vector<size_t> & vec_inliers,
    double & AC_threshold_sq
  );

  static bool estimateRobustPose
  (
    const Frame * frame,
    Hash_Map<MapLandmark* ,IndexT> & vec_putative_matches_3D_frame_idx,
    const size_t min_points,
    Mat3 & R,
    Vec3 & t,
    std::vector<size_t> & vec_inliers,
    double & AC_threshold
  );


  // Find transformation between cameras cam_k and cam_ref: (cam_k)_T_(cam_ref)

  // Look from both directions and find the first intersection
  static bool getRelativeCameraTransformation
  (
    Frame * cam_k,
    Frame * cam_ref,
    Mat4 & T
  );

  // Expess point pt with ref_point_frame in new_point_frame (result in new_pt)
  // For now only with global eculidean representation of points and global represenation of poitns
  static bool getRelativePointPosition(const Vec3 & pt, const Frame * ref_point_frame, Vec3 & pt_new, const Frame * new_point_frame );

  static bool computeFundamentalMatrix(Frame * frame_1, Frame * frame_2, Mat3 & F);

};



}
}

