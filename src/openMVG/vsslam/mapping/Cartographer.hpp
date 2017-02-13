
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
    std::shared_ptr<VSSLAM_Bundle_Adjustment> BA_obj;

    // ---------------------
    // -- Global map
    // ---------------------
    // Keyframes (indexed by frame_id)
    MapFrames keyframes;
    // Landmarks (indexed by landmark_id)
    MapLandmarks structure;
    // Keep track of used landmarkIds
    size_t next_free_landmark_id = 0;

    MAP_POINT_TYPE map_point_type = MAP_POINT_TYPE::EUCLIDEAN;
    MAP_CAMERA_TYPE map_camera_type = MAP_CAMERA_TYPE::ABSOLUTE;

    bool map_initialized = false;

    // ---------------------
    // -- Local map
    // ---------------------
    // Keyframes (indexed by frame_id)
    MapFrames tmp_keyframes;
    std::deque<MapLandmark> tmp_structure;

    // Camera intrinsics  (indexed by camera_id)
    Intrinsics cam_intrinsics;

    // Feature extractor used
    Abstract_FeatureExtractor * feature_extractor_ = nullptr;;

  public:
    //VSSLAM_Data slam_data;

    Cartographer();

    void setFeatureExtractor(Abstract_FeatureExtractor * f_extractor)
     {
       feature_extractor_ = f_extractor;
     }
    // ------------------------------
    // -- Get properties of map and frame representation
    // ------------------------------
    MAP_POINT_TYPE & getMapPointType()
    {
      return map_point_type;
    }
    MAP_CAMERA_TYPE & getCameraPointType()
    {
      return map_camera_type;
    }

    // ------------------------------
    // -- Local Map manipulation
    // ------------------------------
    bool isMapInitialized(){
      return map_initialized;
    }
    void addLocalKeyframe(const std::shared_ptr<Frame> & frame);
    void addLocalMapPoints();
    // ------------------------------
    // -- Map manipulation
    // ------------------------------
    size_t getNextFreeLandmarkId()
    {
      return next_free_landmark_id++;
    }

    void initializeMap
    (
      std::shared_ptr<Frame> & frame_1,
      std::shared_ptr<Frame> & frame_2,
      std::vector<std::pair<Vec3, std::deque<std::pair<Frame*,size_t> > > > & vec_new_pts_3D_obs
    );

    void clearMap();

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
    // -- Add data to map
    // ------------------------------

    void addCameraData(const size_t & cam_id, IntrinsicBase * cam_intrinsic);
    void addKeyframeToMap(const std::shared_ptr<Frame> & frame);
    MapLandmark * newLandmarkToMap();
    void addObservationsToLandmarks(Frame * frame);

    void AddLandmarksToMap
    (
      std::vector<std::pair<Vec3, std::deque<std::pair<Frame*,size_t> > > > & vec_new_pts_3D_obs
    );

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




};

}
}
