// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_BA_CERES_DRIFT_FUNCTOR_HPP
#define OPENMVG_SFM_DATA_BA_CERES_DRIFT_FUNCTOR_HPP

#include "openMVG/cameras/cameras.hpp"
#include "ceres/ceres.h"
#include "ceres/rotation.h"

//--
//- Define ceres Cost_functor for each OpenMVG camera model
//--

namespace openMVG {
namespace sfm {


/// Decorator used to Weight a given cost camera functor
/// i.e useful to weight GCP (Ground Control Points)
template <typename CostFunctor>
struct WeightedDriftCostFunction
{
  WeightedDriftCostFunction() :weight_(1.0) {}

  explicit WeightedDriftCostFunction
  (
    CostFunctor * func,
    const double weight
  )
    :functor_(func), weight_(weight)
  {}

  template <typename T>
  bool operator()
  (
    const T* const pos_3dpoint_A,
    const T* const pos_3dpoint_B,
    T* out_residuals
  ) const
  {
    if (functor_->operator()(pos_3dpoint_A,pos_3dpoint_B, out_residuals))
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
 *  Data parameter blocks are the following <2,3,6,3>
 *  - 2 => dimension of the residuals,
 *  - 3 => the intrinsic data block [focal, principal point x, principal point y],
 *  - 6 => the camera extrinsic data block (camera orientation and position) [R;t],
 *         - rotation(angle axis), and translation [rX,rY,rZ,tx,ty,tz].
 *  - 3 => a 3D point data block.
 *
 */
struct ResidualErrorFunctor_Drift_Point
{
  ResidualErrorFunctor_Drift_Point(const double* const norm_scale)
  :m_pointcloud_norm_scale(norm_scale)
  {
  }

  /**
   * @param[in] pos_3dpoint_A
   * @param[in] pos_3dpoint_B
   * @param[out] out_residuals
   */
  template <typename T>
  bool operator()(
    const T* const pos_3dpoint_A,
    const T* const pos_3dpoint_B,
    T* out_residuals) const
  {

    // Apply focal length and principal point to get the final image coordinates
    const T error_3d_x = (pos_3dpoint_A[0] - pos_3dpoint_B[0])*(pos_3dpoint_A[0] - pos_3dpoint_B[0]);
    const T error_3d_y = (pos_3dpoint_A[1] - pos_3dpoint_B[1])*(pos_3dpoint_A[1] - pos_3dpoint_B[1]);
    const T error_3d_z = (pos_3dpoint_A[2] - pos_3dpoint_B[2])*(pos_3dpoint_A[2] - pos_3dpoint_B[2]);
    

    // Compute and return the error is the difference between the predicted
    //  and observed position
    out_residuals[0] = error_3d_x/(*m_pointcloud_norm_scale);
    out_residuals[1] = error_3d_y/(*m_pointcloud_norm_scale);
    out_residuals[2] = error_3d_z/(*m_pointcloud_norm_scale);
    
    return true;
  }

  static int num_residuals() { return 3; }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction* Create
  (
    const double norm_scale,
    const double weight = 0.0
  )
  {
    if (weight == 0.0)
    {
      return
        (new ceres::AutoDiffCostFunction
          <ResidualErrorFunctor_Drift_Point, 3, 3, 3>(
            new ResidualErrorFunctor_Drift_Point(&norm_scale)));
    }
    else
    {
      return
        (new ceres::AutoDiffCostFunction
          <WeightedDriftCostFunction<ResidualErrorFunctor_Drift_Point>, 3, 3, 3>
          (new WeightedDriftCostFunction<ResidualErrorFunctor_Drift_Point>
            (new ResidualErrorFunctor_Drift_Point(&norm_scale), weight)));
    }
  }



  const double * m_pointcloud_norm_scale; // The 2D observation
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_BA_CERES_DRIFT_FUNCTOR_HPP
