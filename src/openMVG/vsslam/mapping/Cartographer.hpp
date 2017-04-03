
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

//#include <tuple>
/*
#include <string.h>
#include <stdio.h>

#include "slam/LinearSolver_UberBlock.h"
#include "slam/ConfigSolvers.h" // nonlinear graph solvers
#include "slam/3DSolverBase.h" // want C3DJacobians::Quat_to_AxisAngle() and C3DJacobians::AxisAngle_to_Quat()
#include "slam/Sim3_Types.h"
#include "slam/NonlinearSolver_Lambda_DL.h"
#include "slam/Eigenvalues.h"
#include "slam/Sim3SolverBase.h" // C3DJacobians, CBAJacobians, generally useful functions for BA and SE(3), does not need to be included
*/
#include <openMVG/features/features.hpp>
#include <openMVG/numeric/numeric.h>

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"

#include <openMVG/types.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA_slampp.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA_ceres.hpp>

#include <openMVG/vsslam/tracking/PoseEstimation.hpp>
#include <openMVG/vsslam/detection/Abstract_FeatureExtractor.hpp>
#include <deque>

#include <openMVG/sfm/sfm.hpp>

#include <ceres/types.h>
#include <ceres/cost_function.h>
#include <openMVG/vsslam/Frame.hpp>


using namespace openMVG;

namespace openMVG  {
namespace VSSLAM  {


class Cartographer
{
  private:

    // Optimization object
    std::unique_ptr<VSSLAM_Bundle_Adjustment> BA_obj;
    std::unique_ptr<VSSLAM_Bundle_Adjustment> local_BA_obj;

    // Map parameters
    size_t min_obs_per_landmark = 3;
    size_t min_landmark_per_frame = 20;
    size_t init_min_map_pts = 30;
    size_t max_frames_inactive_local_landmark = 5;

    // ---------------------
    // -- Initialization
    // ---------------------
    std::vector<std::shared_ptr<Frame> > init_map_frames;

    // ---------------------
    // -- Global map
    // ---------------------
    bool map_initialized = false;
    // Map settings
    MAP_POINT_TYPE map_point_type = MAP_POINT_TYPE::GLOBAL_EUCLIDEAN;
    MAP_CAMERA_TYPE map_camera_type = MAP_CAMERA_TYPE::GLOBAL;

    // Keyframes (indexed by frame_id)
    MapFrames keyframes;
    // Landmarks (indexed by landmark_id)
    MapLandmarks structure;
    // Camera intrinsics  (indexed by camera_id)
    Intrinsics cam_intrinsics;

    // Keep track of used landmarkIds
    IndexT next_free_landmark_id = 0;

    // ---------------------
    // -- Local map
    // ---------------------
    Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> > tmp_structure;

    IndexT local_p_id = 0;
    // Feature extractor used
    Abstract_FeatureExtractor * feature_extractor_ = nullptr;;

  public:
    // Counter of insertions
    size_t step_id=0;
    Cartographer();


    void setFeatureExtractor(Abstract_FeatureExtractor * f_extractor)
     {
       feature_extractor_ = f_extractor;
     }
    // ------------------------------
    // -- Get properties of map and frame representation
    // ------------------------------
    const MAP_POINT_TYPE & getMapPointType() const
    {
      return map_point_type;
    }
    const MAP_CAMERA_TYPE & getCameraPointType() const
    {
      return map_camera_type;
    }

    bool isMapInitialized()
    {
      return map_initialized;
    }

    void setCeresLocalBA()
    {
      std::cout<<"Set local BA: ceres\n";
      VSSLAM_Bundle_Adjustment_Ceres::BA_options_Ceres options;
      options.linear_solver_type_ = ceres::DENSE_SCHUR;
      local_BA_obj = std::unique_ptr<VSSLAM_Bundle_Adjustment>(new VSSLAM_Bundle_Adjustment_Ceres(options));
    }

    void setCeresGlobalBA()
    {
      VSSLAM_Bundle_Adjustment_Ceres::BA_options_Ceres options;
      options.linear_solver_type_ = ceres::DENSE_SCHUR;
      BA_obj = std::unique_ptr<VSSLAM_Bundle_Adjustment>(new VSSLAM_Bundle_Adjustment_Ceres(options));
    }

    void setSlamPPLocalBA()
    {
      VSSLAM_Bundle_Adjustment_SlamPP::BA_options_SlamPP options;
      local_BA_obj = std::unique_ptr<VSSLAM_Bundle_Adjustment>(new VSSLAM_Bundle_Adjustment_SlamPP(options));
    }

    void setSlamPPGlobalBA()
    {
      VSSLAM_Bundle_Adjustment_SlamPP::BA_options_SlamPP options;
      BA_obj = std::unique_ptr<VSSLAM_Bundle_Adjustment>(new VSSLAM_Bundle_Adjustment_SlamPP(options));
    }

    bool optimizePose
    (
      Frame * frame_i,
      Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_i_idx,
      bool b_use_loss_function = true
    )
    {
      return local_BA_obj->OptimizePose(frame_i, matches_3D_pts_frame_i_idx,b_use_loss_function);
    }
    bool optimizePose
    (
      Frame * frame_i,
      bool b_use_loss_function = true
    )
    {
      return local_BA_obj->OptimizePose(frame_i,b_use_loss_function);
    }
    bool optimizeLocal
    (
      Frame * frame_i,
      std::vector<std::unique_ptr<MapLandmark> > & vec_triangulated_pts,
      bool b_use_loss_function = true
    )
    {
      return local_BA_obj->OptimizeLocal(frame_i, vec_triangulated_pts,b_use_loss_function);
    }


    // ------------------------------
    // -- Set basic map parameters
    // ------------------------------
    void setMapParameters
    (
      const size_t & min_obs_per_l = 3,
      const size_t & min_landmark_per_f = 3,
      const size_t & init_min_map_p = 30,
      size_t max_frames_inactive_local_l = 5
    )
    {
      min_obs_per_landmark = min_obs_per_l;
      min_landmark_per_frame = min_landmark_per_f;
      init_min_map_pts = init_min_map_p;
      max_frames_inactive_local_landmark = max_frames_inactive_local_l;
    }

    // ------------------------------
    // -- Map Initialization
    // ------------------------------
    // Return true if global map is initialized
    const bool & isMapInitialized() const
    {
      return map_initialized;
    }
    // Clear any data needed only for the initialization of global map
    void clearInitializationData();
    // Add Frames/Landmarks/Observations in the initialization stage
    // Once we have enough frames we try to initialize map
    // Return false if initialization failed (not enough points can be added to the system)
    // In false tracking should start from the beginning
    bool initializationAddStep
    (
      std::shared_ptr<Frame> & frame,
      std::vector<std::unique_ptr<MapLandmark> > * vec_new_pts_3D_obs
    );




    // Determine what is the min number of observations landmark need to be considered for global map
    // Observations can be from any frame (not only keyframes)
    double findMinLandmarkQualityForGlobalMapInitialization();
    double findMinLandmarkQualityForDefinedInGlobalMap
    (
      Frame * frame,
      std::vector<std::unique_ptr<MapLandmark> > * vec_new_pts_3D_obs
    );



    // ------------------------------
    // -- Map manipulation
    // ------------------------------
    void clearAllMapData();

    // Return false if we can not add frame to the map -> something went wrong
    bool addStep
    (
      std::shared_ptr<Frame> & frame,
      std::vector<std::unique_ptr<MapLandmark> > * vec_new_pts_3D_obs
    );

    // Check if enough points seen by frame are in global map or in local and will become global
    // by introducing this frame
    bool isFrameDefinedInGlobalMap(Frame * frame)
    {
      if (!isMapInitialized())
        return false;

      size_t n_defined_pts = 0;
      for (MapLandmark * lm : frame->map_points_)
      {
        // we use min_obs -1 because this frame will ad and observation as well
        if (lm && (lm->isActive() || lm->isValidByConnectivityDegree(min_obs_per_landmark-1)))
        {
          n_defined_pts++;
          if (n_defined_pts == min_landmark_per_frame)
            return true;
        }
      }
      return false;
    }



    // ------------------------------
    // -- Frame manipulation
    // ------------------------------
    void addFrameToGlobalMap(const std::shared_ptr<Frame> & frame, bool b_frame_fixed);


    // ------------------------------
    // -- Landmark manipulation
    // ------------------------------
    IndexT getNextFreeLandmarkId()
    {
      return next_free_landmark_id++;
    }
    void updateBestMapPointDescriptor(MapLandmark * ml);
    void addLandmarksToStructure
    (
      Frame * frame,
      std::vector<std::unique_ptr<MapLandmark> > & new_pts,
      const size_t & min_degree_connectivity
    );
    MapLandmark * addLandmarkToLocalMap(std::unique_ptr<MapLandmark> & lm);
    MapLandmark * addLandmarkToGlobalMap(std::unique_ptr<MapLandmark> & lm);

    void eliminateInactiveLocalLandmarks();

    // ------------------------------
    // -- Observation manipulation
    // ------------------------------
    void addObservationsToLandmarks
    (
      Frame * frame,
      const size_t & min_degree_connectivity
    );



    // ------------------------------
    // -- Camera/Landmark representations
    // ------------------------------

    Mat34 getCameraProjectionMatrix(Frame * frame, Frame * frame_ref);

    void verifyLocalLandmarks(Frame *);
    void decreaseLocalFrameCount(Frame * frame);
    void increaseLocalFrameCount(Frame * frame);


    bool exportSceneToPly
    (
      const std::string & filename,
      sfm::ESfM_Data flags_part
    );




    // ------------------------------
    // -- Map manipulation
    // ------------------------------

    Frame * getKeyframeAtPos(size_t i)
    {
      MapFrames::iterator it = keyframes.begin();
      std::advance(it,i);
      return it->second.get();
    };

    size_t getNumberOfKeyframes()
    {
      return keyframes.size();
    };



    // ------------------------------
    // -- Local map
    // ------------------------------
    void getLocalMapPoints
    (
      Frame * currentFrame,
      std::vector<Frame*> & local_frames,
      std::vector<MapLandmark*> & local_points
    );

    void resetFlagLocalMapPoints(std::vector<MapLandmark*> & local_points);
    // void updateLocalMap(Frame * currentFrame);


    void addObservationToIncSystem(MapLandmark * map_point, MapObservation * map_observation)
    {
      BA_obj->addObservationToGlobalSystem(map_point,map_observation);
    };
    void addLandmarkToIncSystem(MapLandmark * map_point)
    {
      BA_obj->addLandmarkToGlobalSysyem(map_point);
    };
    void addFrameToIncSystem(Frame * frame, bool b_frame_fixed = false)
    {
      BA_obj->addFrameToGlobalSystem(frame,b_frame_fixed);
    };
    bool optimizeIncSystem()
    {
      bool b_ok = BA_obj->optimizeGlobal(keyframes,structure);
      std::cout<<"Cartographer: [BA] Global optimization step: "<<step_id<<" was: "<<b_ok<<"!\n";
      return b_ok;
    };


};

}
}
