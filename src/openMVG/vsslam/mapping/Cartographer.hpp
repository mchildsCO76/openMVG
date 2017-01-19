
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once


#include <openMVG/features/features.hpp>
#include <openMVG/numeric/numeric.h>

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"

#include <openMVG/types.hpp>
#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA.hpp>
#include <openMVG/vsslam/optimization/VSSLAM_data_BA_ceres.hpp>
#include <deque>


#include <ceres/types.h>
#include <ceres/cost_function.h>

using namespace openMVG;

namespace openMVG  {
namespace VSSLAM  {

class Cartographer
{
  private:
    VSSLAM_Data slam_data;
    std::shared_ptr<VSSLAM_Bundle_Adjustment> BA_obj;
  public:
    Cartographer();

    void addCameraIntrinsicData(const size_t & cam_id, IntrinsicBase * cam_intrinsic);
    void addObservationToLandmark(const size_t & landmark_id, const std::shared_ptr<Frame> & frame, const size_t featId);
    void addKeyFrameToMap(const std::shared_ptr<Frame> & frame);
    void addLandmarkToMap(MapLandmark & landmark);

    bool initializeMap
    (
      std::shared_ptr<Frame> & frame_1,
      std::shared_ptr<Frame> & frame_2,
      Hash_Map<size_t,size_t> & matches_2_1_idx,
      std::vector<size_t> & inliers
    );

    void clearMap();
  };




  /*


    bool initialization
    (
        std::shared_ptr<Frame> & frame_1,
        std::shared_ptr<Frame> & frame_2,
        Hash_Map<size_t,size_t> & matches_2_1_idx,
        std::vector<size_t> & inliers
    )
    {
      std::cout<<"MAP: Initialization\n";

      // Set owner of second cam if we have relative camera system
      if (mapCameraType == MapLandmarks::MAP_CAMERA_TYPE::RELATIVE)
      {
          frame_2->setReferenceFrame(frame_1);
      }

      // -------------------
      // -- Do tiny BA to refine initial poses
      // -------------------

      Mat34 & P1 = frame_1->getCameraMatrix();
      Mat34 & P2 = frame_2->getCameraMatrix();

      std::vector<Vec3> system_points;
      // Add points to the system
      for (size_t k : inliers)
      {
        Hash_Map<size_t,size_t>::iterator m_iter = matches_2_1_idx.begin();
        std::advance(m_iter,k);

        // 3D point
        Vec3 X;
        // Position of detected features
        const Vec2
          & x1_ = frame_1->getFeaturePosition(m_iter->second),
          & x2_ = frame_2->getFeaturePosition(m_iter->first);
        // Triangulate results
        TriangulateDLT(P1, x1_, P2, x2_, &X);
        system_points.push_back(X);
      }

      //BA

      // -------------------
      // -- Go through results and add poses and cameras
      // -------------------
      // Reset map
      resetMap();
      // Update camera poses

      // Create two keyframes
      addKeyFrameToMap(frame_1);
      addKeyFrameToMap(frame_2);



      // Go through resulting points and see if they are inliers
      size_t pt_i = 0;
      for (size_t k : inliers)
      {
        Vec3 &X = system_points[pt_i];

        if (frame_1->pose_.depth(X) > 0)


      // Add all inliers to the map

        // Normals from each camera center
        Vec3 n_1;
        Vec3 n_2;


        // Get position of both features
        const Vec2
          & x1_ = frame_1->getFeaturePosition(m_iter->second),
          & x2_ = frame_2->getFeaturePosition(m_iter->first);

        TriangulateDLT(P1, x1_, P2, x2_, &X);


          //Vec3 O2 = -R2.transpose() * t2; // Origin of Cam 2


          //float cosParallax;

          for (size_t k = 0; k < vec_inliers.size(); ++k)
          {
            const Vec2
              & x1_ = x1.col(vec_inliers[k]),
              & x2_ = x2.col(vec_inliers[k]);
            // Test if point is front to the two cameras.
            if (Depth(R1, t1, X) < 0 || Depth(R2, t2, X) < 0)
            {
              continue;
            }

            // Check if the parallax between points is sufficient
            n_1 = (X - O1);
            n_2 = (X - O2);
            cosParallax = n_1.dot(n_2)/(n_1.norm() * n_2.norm());

            if (cosParallax < 0.99998)
            {
              continue;
            }

            // Check reprojection error at frame_1
            if ((x1_ - Vec3(P1*X.homogeneous()).hnormalized()).squaredNorm() > E_thresh )
              continue;

            // Check reprojection error at frame_2
            if ((x2_ - Vec3(P2*X.homogeneous()).hnormalized()).squaredNorm() > E_thresh )
              continue;

            ++f[i];
          }




        // Determine point position
        switch (mapPointType)
        {
          case MAP_POINT_TYPE::EUCLIDEAN:
            // Points in world coordinate system
            l.pt_ = *pt3D_iter;
            // Add observations to point
            l.addObservationToLandmark(frame_1,match_2_1_idx.second);
            l.addObservationToLandmark(frame_2,match_2_1_idx.first);
          break;
          case MAP_POINT_TYPE::INV_DIST:
            // Transform point from frame_1 to frame_2
            if (!getRelativePointPosition(*pt3D_iter, frame_1,l.pt_,frame_2))
              continue;
            // Add observations to point
            l.addObservationToLandmark(frame_1,match_2_1_idx.second);
            l.addObservationToLandmark(frame_2,match_2_1_idx.first);
          break;
        }
        // Add landmark to the map
        addLandmarkToMap(l);
      }


      //Create MapPoint.
      //pMP->ComputeDistinctiveDescriptors();
      //pMP->UpdateNormalAndDepth();
*/

}
}
