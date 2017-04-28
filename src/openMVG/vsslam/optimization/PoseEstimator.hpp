// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "openMVG/numeric/numeric.h"
#include "openMVG/cameras/cameras.hpp"
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
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"

#include <openMVG/vsslam/vsslam_parameters.hpp>
#include <openMVG/vsslam/vsslam_data.hpp>
#include <openMVG/vsslam/system/Frame.hpp>
#include <openMVG/vsslam/system/Camera.hpp>

#include <openMVG/vsslam/display/vsslam_time_stats.hpp>


extern openMVG::vsslam::VSSLAM_Time_Stats time_data;

using namespace openMVG;
using namespace openMVG::robust;
using namespace openMVG::geometry;

namespace openMVG {
namespace vsslam {


class Camera;
class Frame;

class PoseEstimator
{
public:
  static const size_t ACRANSAC_ITER = 2048;

  static double computeCosParallaxBetweenRays
  (
    const Vec3 & ray_1,
    const Vec3 & ray_2
  );

  static bool computeH
  (
    const Mat & x1,
    const Mat & x2,
    const std::pair<size_t,size_t> & size_img_1,
    const std::pair<size_t,size_t> & size_img_2,
    const double & f_max_thresh,
    Mat3 & H_21,  //x2 = H21 * x1
    std::pair<double, double> & ACRansacInfo,
    std::vector<unsigned int> & vec_inliers
  );

  static bool computeE
  (
    const Mat3 & K1,
    const Mat3 & K2,
    const Mat & x1,
    const Mat & x2,
    const std::pair<size_t,size_t> & size_img_1,
    const std::pair<size_t,size_t> & size_img_2,
    const double & f_max_thresh,
    Mat3 & E_21,
    std::pair<double, double> & ACRansacInfo,
    std::vector<unsigned int> & vec_inliers
  );

  static bool estimateRobustRelativePose_HE_Pinhole
  (
    const Frame * frame_1,
    const Frame * frame_2,
    const openMVG::matching::IndMatches & vec_putative_matches_1_2_idx,
    Mat4 & T,
    std::vector<uint32_t> & vec_inliers,
    double & f_model_thresh,
    VSSLAM_Parameters * params
  );

  static bool estimateRobustPose_Pinhole
  (
    const Frame * frame,
    Hash_Map<MapLandmark* ,IndexT> & matches_landmarks_frame,
    Mat4 & T,
    std::vector<uint32_t> & vec_inliers,
    double & f_model_thresh,
    VSSLAM_Parameters * param
  );

  static bool getRelativePointPosition
  (
    const Vec3 & pt_3D_frame,
    const Frame * frame_reference,
    Vec3 & pt_3D_frame_new,
    const Frame * frame_new_reference
  );

  static bool computeFundamentalMatrix
  (
    const Frame * frame_1,
    const Frame * frame_2,
    Mat3 & F_21
  );

  static bool getRelativeCameraTransformation
  (
    const Frame * frame,
    const Frame * frame_reference,
    Mat4 & T
  );
};

}
}
