
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/sfm/sfm_data_BA.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA.hpp>
#include <openMVG/numeric/numeric.h>
#include <ceres/ceres.h>
#include <ceres/types.h>
#include <ceres/cost_function.h>
#include <openMVG/vsslam/Frame.hpp>

using namespace openMVG;

namespace openMVG {
namespace VSSLAM {

ceres::CostFunction * IntrinsicsToCostFunction
(
  IntrinsicBase * intrinsic,
  const Vec2 & observation,
  const Eigen::Matrix<double, 2, 2> & inf_matrix,
  const double weight = 0.0
);




class VSSLAM_Bundle_Adjustment_Ceres : public VSSLAM_Bundle_Adjustment
{
public:
  struct BA_options_Ceres : public BA_options
  {
    unsigned int nb_threads_;
    bool b_ceres_summary_;
    ceres::LinearSolverType linear_solver_type_;
    ceres::PreconditionerType preconditioner_type_;
    ceres::SparseLinearAlgebraLibraryType sparse_linear_algebra_library_type_;
    double parameter_tolerance_;
    bool b_use_loss_function_;

    BA_options_Ceres
    (
      const MAP_CAMERA_TYPE map_camera_type = MAP_CAMERA_TYPE::GLOBAL,
      const MAP_POINT_TYPE map_landmark_type = MAP_POINT_TYPE::GLOBAL_EUCLIDEAN,
      const bool b_verbose = true, bool bmultithreaded = true
    ) : BA_options(map_camera_type, map_landmark_type, b_verbose),
    nb_threads_(1),
    parameter_tolerance_(1e-8), //~= numeric_limits<float>::epsilon()
    b_use_loss_function_(true)
    {
      #ifdef OPENMVG_USE_OPENMP
        nb_threads_ = omp_get_max_threads();
      #endif // OPENMVG_USE_OPENMP
      if (!bmultithreaded)
        nb_threads_ = 1;

      b_ceres_summary_ = false;

      // Default configuration use a DENSE representation
      linear_solver_type_ = ceres::DENSE_SCHUR;
      preconditioner_type_ = ceres::JACOBI;
      // If Sparse linear solver are available
      // Descending priority order by efficiency (SUITE_SPARSE > CX_SPARSE > EIGEN_SPARSE)
      if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
      {
        sparse_linear_algebra_library_type_ = ceres::SUITE_SPARSE;
        linear_solver_type_ = ceres::SPARSE_SCHUR;
      }
      else
      {
        if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
        {
          sparse_linear_algebra_library_type_ = ceres::CX_SPARSE;
          linear_solver_type_ = ceres::SPARSE_SCHUR;
        }
        else
        if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
        {
          sparse_linear_algebra_library_type_ = ceres::EIGEN_SPARSE;
          linear_solver_type_ = ceres::SPARSE_SCHUR;
        }
      }
    }
  };
private:
  BA_options_Ceres options_;
  ceres::Problem problem_;
  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics_;
  Hash_Map<IndexT, std::vector<double> > map_poses_;

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction_;



public:
  VSSLAM_Bundle_Adjustment_Ceres
  (
    BA_options_Ceres options = BA_options_Ceres()
  );

  BA_options_Ceres & getOptions();


  bool OptimizePose
  (
    Frame * frame_i,
    Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_i_idx,
    bool b_use_loss_function = true
  )override;

  bool OptimizePose
  (
    Frame * frame_i,
    bool b_use_loss_function = true
  )override;

  bool OptimizeLocal
  (
    Frame * frame_i,
    std::vector<std::unique_ptr<MapLandmark> > & vec_triangulated_pts,
    bool b_use_loss_function = true
  )override;


  bool addObservationToGlobalSystem(MapLandmark * map_point, MapObservation * map_observation)override;
  bool addLandmarkToGlobalSysyem(MapLandmark * map_point)override;
  bool addFrameToGlobalSystem(Frame * frame, bool b_frame_fixed = false)override;
  bool optimizeGlobal(MapFrames & map_frames, MapLandmarks & map_landmarks)override;


};

}
}
