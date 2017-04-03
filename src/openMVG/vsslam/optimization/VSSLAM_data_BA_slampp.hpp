
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/sfm/sfm_data_BA.hpp>
#include <openMVG/numeric/numeric.h>

#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA.hpp>


#include <openMVG/vsslam/optimization/slampp_optimizers/SlamPP_Optimizer_Sim3_gXYZ_gXYZ.hpp>


using namespace openMVG;

namespace openMVG {
namespace VSSLAM {


class VSSLAM_Bundle_Adjustment_SlamPP : public VSSLAM_Bundle_Adjustment
{
  public:
    struct BA_options_SlamPP : public BA_options
    {
      size_t undefined_cam_id = std::numeric_limits<size_t>::max() / 2 + 1;

      bool b_use_schur_ = true;
      bool b_do_marginals_ = false;
      bool b_do_icra_style_marginals_ = false;
      size_t n_skip_optimization = 0; //Optimize every n included cameras (0 - means every)

      /*
       The default nonlinear solve threshold is 0.005, the default final number of iterations
       is 5, the default update threshold is 10^-4 and the default trust radius is 0.2.
       */

      int n_max_inc_iters = 1;  //max-nonlinear-solve-iters
      int n_max_final_iters = 5;  //max-final-nonlinear-solve-iters
      double f_inc_nlsolve_thresh = .005; //nonlinear-solve-error-thresh
      double f_final_nlsolve_thresh = .005; //final-nonlinear-solve-error-thresh
      double f_trust_radius = 0.2;  // trust-radius
      bool b_trust_radius_persistent = false;
      double f_update_thresh = 0.0001;  //dx-threshold  (using update threshold 0, disabling update threshold)

      bool b_all_batch = false;

      BA_options_SlamPP
      (
        const MAP_CAMERA_TYPE map_camera_type = MAP_CAMERA_TYPE::GLOBAL,
        const MAP_POINT_TYPE map_landmark_type = MAP_POINT_TYPE::GLOBAL_EUCLIDEAN,
        bool b_verbose = false,
        bool b_use_schur = true,
        bool b_do_marginals = false,
        bool b_do_icra_style_marginals = false
      ) : BA_options(map_camera_type, map_landmark_type, b_verbose)
      {
        b_use_schur_ = b_use_schur;
        b_do_marginals_ = b_do_marginals; //~= numeric_limits<float>::epsilon()
        b_do_icra_style_marginals_ = b_do_icra_style_marginals;


      }
    };
  private:

    BA_options_SlamPP options_;
    std::unique_ptr<SlamPP_Optimizer> problem_;
    Hash_Map<IndexT, std::pair<size_t,double *> > map_poses_;
    Hash_Map<IndexT, std::pair<size_t,double *> > map_landmarks_;

    size_t next_vertex_idx_slampp_ = 0;
    size_t next_const_vertex_idx_slampp_ = std::numeric_limits<size_t>::max(); // Constant vertices have decreasing ids starting with max

    size_t getNextVertexId()
    {
      return next_vertex_idx_slampp_++;
    }
    size_t getNextConstVertexId()
    {
      return next_const_vertex_idx_slampp_--;
    }

  public:

    VSSLAM_Bundle_Adjustment_SlamPP
    (
        BA_options_SlamPP options = BA_options_SlamPP()
    );

    BA_options_SlamPP & getOptions();


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
