
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <tuple>
#include <openMVG/features/features.hpp>
#include <openMVG/numeric/numeric.h>

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"

#include <openMVG/types.hpp>
#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA_ceres.hpp>
#include <openMVG/vsslam/tracking/PoseEstimation.hpp>
#include <openMVG/vsslam/tracking/Abstract_FeatureExtractor.hpp>
#include <deque>


#include <ceres/types.h>
#include <ceres/cost_function.h>

using namespace openMVG;

namespace openMVG  {
namespace VSSLAM  {


class Cartographer
{
  private:
    // Counter of insertions
    size_t step_id=0;

    // Optimization object
    std::shared_ptr<VSSLAM_Bundle_Adjustment> BA_obj;

    // Map parameters
    size_t min_obs_per_landmark = 3;
    size_t min_landmark_per_frame = 3;
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
    size_t next_free_landmark_id = 0;

    // ---------------------
    // -- Local map
    // ---------------------
    Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> > tmp_structure;

    // Feature extractor used
    Abstract_FeatureExtractor * feature_extractor_ = nullptr;;

  public:
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


    // Determine what is the min degree of landmark necessary to be able to initialize the map
    size_t findMinLandmarkDegreeForGlobalMapInitialization();
    size_t findMinLandmarkDegreeForDefinedInGlobalMap
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
    void addFrameToGlobalMap(const std::shared_ptr<Frame> & frame);


    // ------------------------------
    // -- Landmark manipulation
    // ------------------------------
    size_t getNextFreeLandmarkId()
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
    MapLandmark * addLandmarkFromLocalToGlobalMap(std::unique_ptr<MapLandmark> & lm);

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

    // void updateLocalMap(Frame * currentFrame);


    void addObservationToIncSystem()
    {};
    void addLandmarkToIncSystem(){};
    void addFrameToIncSystem(){};
    void optimizeIncSystem(){};

};

}
}
