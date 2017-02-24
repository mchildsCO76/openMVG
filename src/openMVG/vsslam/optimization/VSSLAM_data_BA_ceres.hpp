
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/sfm/sfm_data_BA.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA.hpp>
#include <openMVG/numeric/numeric.h>
#include <ceres/types.h>
#include <ceres/cost_function.h>
#include <openMVG/vsslam/Frame.hpp>

using namespace openMVG;

namespace openMVG {
namespace VSSLAM {

ceres::CostFunction * IntrinsicsToCostFunction
(
  cameras::IntrinsicBase * intrinsic,
  const Vec2 & observation,
  const double weight = 0.0
);

class VSSLAM_Bundle_Adjustment_Ceres : public VSSLAM_Bundle_Adjustment
{
  public:
    struct BA_Ceres_options
    {
      bool bVerbose_;
      unsigned int nb_threads_;
      bool bCeres_summary_;
      ceres::LinearSolverType linear_solver_type_;
      ceres::PreconditionerType preconditioner_type_;
      ceres::SparseLinearAlgebraLibraryType sparse_linear_algebra_library_type_;
      double parameter_tolerance_;
      bool bUse_loss_function_;

      BA_Ceres_options(const bool bVerbose = true, bool bmultithreaded = true);
    };
  private:
    VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options ceres_options_;
  public:

    VSSLAM_Bundle_Adjustment_Ceres
    (
      VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options options
    );

    VSSLAM_Bundle_Adjustment_Ceres::BA_Ceres_options & ceres_options();

    /*
    bool OptimizePose
    (
      std::vector<Frame*> * vec_frames,
      Frame * frame_i,
      Hash_Map<MapLandmark *,IndexT> * matches_3D_pts_frame_i_idx,
      std::vector<std::unique_ptr<MapLandmark> > * vec_triangulated_pts
    ) override;
*/
    bool OptimizePose
    (
      Frame * frame_i,
      Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_i_idx
    )override;

    bool OptimizeLocal
    (
      Hash_Map<Frame*, size_t> & tmp_frames,
      Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> > & tmp_structure,
      Frame * frame_i,
      std::vector<std::unique_ptr<MapLandmark> > & vec_triangulated_pts
    )override;

};

}
}
