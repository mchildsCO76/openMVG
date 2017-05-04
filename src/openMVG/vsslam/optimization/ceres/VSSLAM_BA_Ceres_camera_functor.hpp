// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "openMVG/cameras/cameras.hpp"

namespace openMVG {
namespace vsslam {

/// Decorator used to Weight a given cost camera functor
/// i.e useful to weight GCP (Ground Control Points)
template <typename CostFunctor>
struct WeightedCostFunction
{
  WeightedCostFunction() :weight_(1.0) {}

  explicit WeightedCostFunction
  (
    CostFunctor * func,
    const double weight
  )
    :functor_(func), weight_(weight)
  {}

  template <typename T>
  bool operator()
  (
    const T* const cam_intrinsic,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals
  ) const
  {
    if (functor_->operator()(cam_intrinsic, cam_extrinsics, pos_3dpoint, out_residuals))
    {
      // Reweight the residual values
      for (int i = 0; i < CostFunctor::num_residuals(); ++i)
      {
        out_residuals[i] *= T(weight_);
      }
      return true;
    }
    return false;
  }

  ceres::internal::scoped_ptr<CostFunctor> functor_;
  const double weight_;
};


/**
 * @brief Ceres functor to use a Pinhole_Intrinsic (pinhole camera model K[R[t]) and a 3D point.
 *
 *  Data parameter blocks are the following <2,3,7,3>
 *  - 2 => dimension of the residuals,
 *  - 3 => the intrinsic data block [focal, principal point x, principal point y],
 *  - 7 => the camera extrinsic data block (camera orientation and position) [R;t;s],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz;s].
 *  - 3 => a 3D point data block.
 *
 */
struct ResidualErrorFunctor_Pinhole_Intrinsic_Rts
{
  ResidualErrorFunctor_Pinhole_Intrinsic_Rts(const double* const pos_2dpoint)
  :m_pos_2dpoint(pos_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 7 parameters [R;t;s]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];
    const T * cam_s = &cam_extrinsics[6];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply scale
    /*pos_proj[0] *= *cam_s;
    pos_proj[1] *= *cam_s;
    pos_proj[2] *= *cam_s;
*/
    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_u;
    const T projected_y = principal_point_y + focal * y_u;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = projected_x - T(m_pos_2dpoint[0]);
    out_residuals[1] = projected_y - T(m_pos_2dpoint[1]);

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Pinhole_Intrinsic_Rts, 2, 3, 7, 3>(
            new ResidualErrorFunctor_Pinhole_Intrinsic_Rts(observation.data())));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Rts>, 2, 3, 7, 3>
          (new WeightedCostFunction<ResidualErrorFunctor_Pinhole_Intrinsic_Rts>
            (new ResidualErrorFunctor_Pinhole_Intrinsic_Rts(observation.data()), weight)));
    }
  }

  const double * m_pos_2dpoint; // The 2D observation
};

struct Chi2ErrorFunctor_Pinhole_Intrinsic_Rts
{

  Chi2ErrorFunctor_Pinhole_Intrinsic_Rts(const double* const pos_2dpoint,const double* const inf_2dpoint)
  :m_pos_2dpoint(pos_2dpoint), m_sqrt_inf_mat_2dpoint(inf_2dpoint)
  {
  }

  // Enum to map intrinsics parameters between openMVG & ceres camera data parameter block.
  enum {
    OFFSET_FOCAL_LENGTH = 0,
    OFFSET_PRINCIPAL_POINT_X = 1,
    OFFSET_PRINCIPAL_POINT_Y = 2
  };

  /**
   * @param[in] cam_intrinsics: Camera intrinsics( focal, principal point [x,y] )
   * @param[in] cam_extrinsics: Camera parameterized using one block of 7 parameters [R;t;s]:
   *   - 3 for rotation(angle axis), 3 for translation
   * @param[in] pos_3dpoint
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const cam_intrinsics,
    const T* const cam_extrinsics,
    const T* const pos_3dpoint,
    T* out_residuals) const
  {
    //--
    // Apply external parameters (Pose)
    //--

    const T * cam_R = cam_extrinsics;
    const T * cam_t = &cam_extrinsics[3];
    const T * cam_s = &cam_extrinsics[6];

    T pos_proj[3];
    // Rotate the point according the camera rotation
    ceres::AngleAxisRotatePoint(cam_R, pos_3dpoint, pos_proj);

    // Apply scale
    /*pos_proj[0] *= *cam_s;
    pos_proj[1] *= *cam_s;
    pos_proj[2] *= *cam_s;
*/
    // Apply the camera translation
    pos_proj[0] += cam_t[0];
    pos_proj[1] += cam_t[1];
    pos_proj[2] += cam_t[2];

    // Transform the point from homogeneous to euclidean (undistorted point)
    const T x_u = pos_proj[0] / pos_proj[2];
    const T y_u = pos_proj[1] / pos_proj[2];

    //--
    // Apply intrinsic parameters
    //--

    const T& focal = cam_intrinsics[OFFSET_FOCAL_LENGTH];
    const T& principal_point_x = cam_intrinsics[OFFSET_PRINCIPAL_POINT_X];
    const T& principal_point_y = cam_intrinsics[OFFSET_PRINCIPAL_POINT_Y];

    // Apply focal length and principal point to get the final image coordinates
    const T projected_x = principal_point_x + focal * x_u;
    const T projected_y = principal_point_y + focal * y_u;

    // Compute and return the error is the difference between the predicted
    //  and observed position
    const T e_x = projected_x - m_pos_2dpoint[0];
    const T e_y = projected_y - m_pos_2dpoint[1];

    // Eigen matrix m_inf_2dpoint is stored in ColumnMajor order (default by Eigen)
    out_residuals[0] = T(m_sqrt_inf_mat_2dpoint[0]) * e_x + T(m_sqrt_inf_mat_2dpoint[2]) * e_y;
    out_residuals[1] = T(m_sqrt_inf_mat_2dpoint[1]) * e_x + T(m_sqrt_inf_mat_2dpoint[3]) * e_y;

    return true;
  }

  static int num_residuals() { return 2; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const Vec2 & observation,
    const Eigen::Matrix<double, 2, 2> & sqrt_inf_matrix,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <Chi2ErrorFunctor_Pinhole_Intrinsic_Rts, 2, 3, 7, 3>(
            new Chi2ErrorFunctor_Pinhole_Intrinsic_Rts(observation.data(),sqrt_inf_matrix.data())));

    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedCostFunction<Chi2ErrorFunctor_Pinhole_Intrinsic_Rts>, 2, 3, 7, 3>
          (new WeightedCostFunction<Chi2ErrorFunctor_Pinhole_Intrinsic_Rts>
            (new Chi2ErrorFunctor_Pinhole_Intrinsic_Rts(observation.data(),sqrt_inf_matrix.data()), weight)));

    }
  }

  const double * m_pos_2dpoint; // The 2D observation
  const double * m_sqrt_inf_mat_2dpoint; // The information of 2D observation

};

}
}
