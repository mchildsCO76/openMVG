
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/vsslam/matching/Abstract_FeatureMatcher.hpp>
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"

namespace openMVG  {
namespace VSSLAM  {

class Feat_Matcher_CascadeHashing : public Abstract_FeatureMatcher
{
private:

  void purgeCandidateMatches
  (
    std::vector<IndexT> & vec_putative_matches_1_2_idx,
    matching::IndMatches & vec_matches_1_2_idx
  )
  {
    // Check that two are not matching the same point
    bool bOkMatch = true;
    for (IndexT i=0; i<vec_putative_matches_1_2_idx.size(); ++i)
    {
      if (vec_putative_matches_1_2_idx[i] != UndefinedIndexT)
      {
        for (IndexT j=i+1; j<vec_putative_matches_1_2_idx.size(); ++j)
        {
          // if value is doubled we delete both matches (mathces have to be unique)
          if (vec_putative_matches_1_2_idx[i] == vec_putative_matches_1_2_idx[j])
          {
            bOkMatch = false;
            vec_putative_matches_1_2_idx[j] = UndefinedIndexT;
          }
        }
        if (bOkMatch)
        {
          // Match
          vec_matches_1_2_idx.emplace_back(i,vec_putative_matches_1_2_idx[i]);
        }
        else
        {
          // reset duplicate flag
          bOkMatch = true;

        }
      }
    }
  }

  void purgeCandidateMatches
  (
    const std::vector<MapLandmark *> & vec_3D_pts,
    std::vector<IndexT> & putative_matches_3D_frame_idx,
    Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_idx
  )
  {
    // Check that two are not matching the same point
    bool bOkMatch = true;
    for (size_t i=0; i<putative_matches_3D_frame_idx.size(); ++i)
    {
      if (putative_matches_3D_frame_idx[i] != UndefinedIndexT)
      {
        for (size_t j=i+1; j<putative_matches_3D_frame_idx.size(); ++j)
        {
          // if value is doubled we delete both matches (mathces have to be unique)
          if (putative_matches_3D_frame_idx[i] == putative_matches_3D_frame_idx[j])
          {
            bOkMatch = false;
            putative_matches_3D_frame_idx[j] = UndefinedIndexT;
          }
        }
        if (bOkMatch)
        {
          // Match
          matches_3D_pts_frame_idx[vec_3D_pts[i]] = putative_matches_3D_frame_idx[i];
        }
        else
        {
          // reset duplicate flag
          bOkMatch = true;
        }
      }
    }
  }

  void matchingWithCascadeHashing
  (
    const Frame * frame_1,
    const Frame * frame_2,
    matching::IndMatches & vec_putative_matches,
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  )
  {
    const std::string descriptor_type = frame_1->getRegions()->Type_id();

    if (descriptor_type == typeid(float).name())
    {
      matchingWithCascadeHashing_All_All_2D_2D<float>(frame_1,frame_2,vec_putative_matches,desc_ratio,max_desc_d);
    }
    else if (descriptor_type == typeid(double).name())
    {
      matchingWithCascadeHashing_All_All_2D_2D<double>(frame_1,frame_2,vec_putative_matches,desc_ratio,max_desc_d);
    }
    else if (descriptor_type == typeid(unsigned char).name())
    {
      matchingWithCascadeHashing_All_All_2D_2D<unsigned char>(frame_1,frame_2,vec_putative_matches,desc_ratio,max_desc_d);
    }
  }

  template <typename ScalarT>
  void matchingWithCascadeHashing_All_All_2D_2D
  (
    const Frame * frame_1,
    const Frame * frame_2,
    matching::IndMatches & vec_putative_matches,
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  )
  {
    // Compute the zero mean descriptor that will be used for hashing (one for all the image regions)
    features::Regions * const frame_1_regions = frame_1->getRegions();
    features::Regions * const frame_2_regions = frame_2->getRegions();

    const size_t dimension = frame_1_regions->DescriptorLength();
    Eigen::MatrixXf matForZeroMean(2,dimension);

    const ScalarT * tab_1 = reinterpret_cast<const ScalarT *>(frame_1->getRegions()->DescriptorRawData());
    Eigen::Map<Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > mat_1( (ScalarT*)tab_1, frame_1->getNumberOfFeatures(), dimension);
    matForZeroMean.row(0) = matching::CascadeHasher::GetZeroMeanDescriptor(mat_1);

    const ScalarT * tab_2 = reinterpret_cast<const ScalarT *>(frame_2->getRegions()->DescriptorRawData());
    Eigen::Map<Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > mat_2( (ScalarT*)tab_2, frame_2->getNumberOfFeatures(), dimension);
    matForZeroMean.row(1) = matching::CascadeHasher::GetZeroMeanDescriptor(mat_2);


    Eigen::VectorXf zero_mean_descriptor = matching::CascadeHasher::GetZeroMeanDescriptor(matForZeroMean);

    // Compute hashed descriptors (with zero mean vector)
    matching::HashedDescriptions hashed_descriptors_1 = cascade_hasher_.CreateHashedDescriptions(mat_1,zero_mean_descriptor);
    matching::HashedDescriptions hashed_descriptors_2 = cascade_hasher_.CreateHashedDescriptions(mat_2,zero_mean_descriptor);


    matching::IndMatches pvec_indices;
    using ResultType = typename Accumulator<ScalarT>::Type;
    std::vector<ResultType> pvec_distances;
    pvec_distances.reserve(frame_2->getNumberOfFeatures() * 2);
    pvec_indices.reserve(frame_2->getNumberOfFeatures() * 2);

    // Match the query descriptors to the database
    cascade_hasher_.Match_HashedDescriptions<Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, ResultType>(
      hashed_descriptors_2, mat_2,
      hashed_descriptors_1, mat_1,
      &pvec_indices, &pvec_distances);

    std::vector<int> vec_nn_ratio_idx;
    // Filter the matches using a distance ratio test:
    //   The probability that a match is correct is determined by taking
    //   the ratio of distance from the closest neighbor to the distance
    //   of the second closest.

    matching::NNdistanceRatio(
      pvec_distances.begin(), // distance start
      pvec_distances.end(),   // distance end
      2, // Number of neighbor in iterator sequence (minimum required 2)
      vec_nn_ratio_idx, // output (indices that respect the distance Ratio)
      desc_ratio);

    vec_putative_matches.reserve(vec_nn_ratio_idx.size());
    for (size_t k=0; k < vec_nn_ratio_idx.size(); ++k)
    {
      const size_t index = vec_nn_ratio_idx[k];
      // Check if distance of best descriptor is below threshold (useful?)
      if (pvec_distances[index*2] > max_desc_d)
        continue;
      vec_putative_matches.emplace_back(pvec_indices[index*2].j_, pvec_indices[index*2].i_);
    }
    // Remove duplicates
    matching::IndMatch::getDeduplicated(vec_putative_matches);

    // Remove matches that have the same (X,Y) coordinates
    /*matching::IndMatchDecorator<float> matchDeduplicator(vec_putative_matches,
        frame_1_regions->GetRegionsPositions(), frame_2_regions->GetRegionsPositions());
    matchDeduplicator.getDeduplicated(vec_putative_matches);
*/

  }


public:
  matching::CascadeHasher cascade_hasher_;

  Feat_Matcher_CascadeHashing
  (
    Abstract_FeatureExtractor * feat_extractor
  )
  {
    cascade_hasher_.Init(feat_extractor->getDescriptorLength());
  }


  // Match two sets of features by:
  //   for each feature in F_1:
  //     compare descriptors only of features in F_2
  //     that are inside a window around position of feature in F_1


  void matching_AllAll_2D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    matching::IndMatches & vec_putative_matches_1_2_idx,
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    matchingWithCascadeHashing(frame_1,frame_2,vec_putative_matches_1_2_idx,desc_ratio,max_desc_d);
  }

  // Matching all available 3D points with the 2D features of frame_2
  void matching_AllAll_3D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const std::vector<MapLandmark *> & vec_3D_pts,
    const Frame * frame,
    Hash_Map<MapLandmark* ,IndexT> & matches_3D_pts_frame_idx,
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float desc_ratio_2 = desc_ratio*desc_ratio;
    const float max_desc_d_2 = max_desc_d*max_desc_d;

    // 3D data
    const size_t vec_3D_pts_size = vec_3D_pts.size();
    // Frame data
    features::Regions * const frame_regions = frame->getRegions();
    const std::vector<Vec2> & frame_pts = frame->pts_undist_;
    const IntrinsicBase * cam_intrinsic = frame->getCameraIntrinsics();

    // -------------------
    // -- Match features of prev_frame(that are already triangulated and in map) with the features of the new frame
    // -------------------
    std::vector<IndexT>matches_3D_frame_idx(vec_3D_pts_size,UndefinedIndexT);

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t p_i=0; p_i < vec_3D_pts_size; ++p_i)
    {
      MapLandmark * pt_3D = vec_3D_pts[p_i];

      // Check if its triangulated
      if (!pt_3D)
        continue;

      // Already matched
      if (matches_3D_pts_frame_idx.find(pt_3D) != matches_3D_pts_frame_idx.end())
        continue;

      // Project the point to frame 2
      Vec3 pt_3D_frame;
      getRelativePointPosition(pt_3D->X_,pt_3D->ref_frame_,pt_3D_frame,frame);
      if (pt_3D_frame(2) < 0)
        continue;

      // Project the point to image plane of frame 2
      const Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(cam_intrinsic->have_disto()?cam_intrinsic->add_disto(pt_3D_frame.hnormalized()):pt_3D_frame.hnormalized());

      // Check if projection is actually in the image borders
      if (!frame->isPointInFrame(pt_3D_frame_projected))
        continue;

      // Raw descriptors
      void * pt_3D_raw_desc = pt_3D->bestDesc_;
      void * candidate_pt_desc_raw;

      //TODO: Get possible candidates through grid

      // Check for the best and second best match
      IndexT best_idx = UndefinedIndexT;
      IndexT second_best_idx = UndefinedIndexT;
      double best_distance_2 = max_desc_d_2;
      double second_best_distance_2 = std::numeric_limits<double>::infinity();
      double distance_2;

      for (size_t c_i=0; c_i < frame_pts.size(); ++c_i)
      {
        // This point in the frame is already taken
        if (frame->map_points_[c_i])
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance_2
        distance_2 = featureExtractor_->SquaredDescriptorDistance(pt_3D_raw_desc,candidate_pt_desc_raw);

        // Save if in best two
        if (distance_2 < best_distance_2)
        {
          second_best_distance_2 = best_distance_2;
          second_best_idx = best_idx;
          best_distance_2 = distance_2;
          best_idx = c_i;
        }
        else if (distance_2 < second_best_distance_2)
        {
          second_best_distance_2 = distance_2;
          second_best_idx = c_i;
        }
      }

      // Detect best match
      if (best_idx != UndefinedIndexT)
      {
        if (second_best_idx != UndefinedIndexT  && (best_distance_2 / second_best_distance_2) < desc_ratio_2)
        {
          // Best is unique enough
          matches_3D_frame_idx[p_i] = best_idx;
        }
        else
        {
          // Best is unique
          matches_3D_frame_idx[p_i] = best_idx;
        }
      }
    }

    // -------------------
    // -- Purge matches
    // -------------------
    purgeCandidateMatches(vec_3D_pts,matches_3D_frame_idx,matches_3D_pts_frame_idx);

  }

  void matching_Projection_3D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame,
    const std::vector<MapLandmark *> & vec_3D_pts,
    Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_idx,
    const size_t max_dist = 15, //20*20
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float desc_ratio_2 = desc_ratio*desc_ratio;
    const float max_desc_d_2 = max_desc_d * max_desc_d;
    const float max_dist_2 = max_dist * max_dist;

    // TODO: adjust max_d_2 from the viewing angle

    // 3D data
    const size_t vec_3D_pts_size = vec_3D_pts.size();
    // Frame data
    features::Regions * const frame_regions = frame->getRegions();
    const std::vector<Vec2> & frame_pts = frame->pts_undist_;
    const IntrinsicBase * cam_intrinsic = frame->getCameraIntrinsics();

    // -------------------
    // -- Match features of prev_frame(that are already triangulated and in map) with the features of the new frame
    // -------------------
    std::vector<IndexT>matches_3D_frame_idx(vec_3D_pts_size,UndefinedIndexT);

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t p_i=0; p_i < vec_3D_pts_size; ++p_i)
    {
      MapLandmark * pt_3D = vec_3D_pts[p_i];

      // Check if its triangulated
      if (!pt_3D)
        continue;

      // Already matched
      if (matches_3D_pts_frame_idx.find(pt_3D) != matches_3D_pts_frame_idx.end())
        continue;

      // Project the point to frame 2
      Vec3 pt_3D_frame;
      getRelativePointPosition(pt_3D->X_,pt_3D->ref_frame_,pt_3D_frame,frame);
      if (pt_3D_frame(2) < 0)
        continue;

      // Project the point to image plane of frame 2
      const Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(cam_intrinsic->have_disto()?cam_intrinsic->add_disto(pt_3D_frame.hnormalized()):pt_3D_frame.hnormalized());

      // Check if projection is actually in the image borders
      if (!frame->isPointInFrame(pt_3D_frame_projected))
        continue;


      // Compute viewing angle between the average angle of the point and
      Vec3 P0 = pt_3D->getWorldPosition() - frame->O_w_;
      Vec3 Pn = pt_3D->getNormal();

      const float cos_view_angle = P0.dot(Pn) / P0.norm();
      // enlarge the area if the viewing angle is bigger
      float factor = radiusByViewingAngle(cos_view_angle);
      float max_dist_2 = (max_dist * factor) * (max_dist * factor);

      // Raw descriptors
      void * pt_3D_raw_desc = pt_3D->bestDesc_;
      void * candidate_pt_desc_raw;

      //TODO: Get possible candidates through grid

      // Loop through all the candidates and find the best and second best
      IndexT best_idx = UndefinedIndexT;
      IndexT second_best_idx = UndefinedIndexT;
      double best_distance_2 = max_desc_d_2;
      double second_best_distance_2 = std::numeric_limits<double>::infinity();
      double distance_2;

      for (size_t c_i=0; c_i < frame_pts.size(); ++c_i)
      {
        // This point in the frame is already taken
        if (frame->map_points_[c_i])
          continue;


        // Check if points are close enough
        if ((frame_pts[c_i] - pt_3D_frame_projected).squaredNorm() > max_dist_2)
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance_2
        distance_2 = featureExtractor_->SquaredDescriptorDistance(pt_3D_raw_desc,candidate_pt_desc_raw);

        // Save if in best two
        if (distance_2 < best_distance_2)
        {
          second_best_distance_2 = best_distance_2;
          second_best_idx = best_idx;
          best_distance_2 = distance_2;
          best_idx = c_i;
        }
        else if (distance_2 < second_best_distance_2)
        {
          second_best_distance_2 = distance_2;
          second_best_idx = c_i;
        }
      }

      // Detect best match
      if (best_idx != UndefinedIndexT)
      {
        if (second_best_idx != UndefinedIndexT && (best_distance_2 / second_best_distance_2) < desc_ratio_2)
        {
          // Best is unique enough
          matches_3D_frame_idx[p_i] = best_idx;
        }
        else
        {
          // Best is unique
          matches_3D_frame_idx[p_i] = best_idx;
        }
      }
    }

    // -------------------
    // -- Purge matches
    // -------------------
    purgeCandidateMatches(vec_3D_pts,matches_3D_frame_idx,matches_3D_pts_frame_idx);

  }

  void matching_EpipolarLine_2D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    const Mat3 & F_21,
    Hash_Map<IndexT,IndexT> & map_matches_1_2_idx,
    const size_t max_epipolar_dist = 16,
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {

    const float max_epipolar_d_2 = max_epipolar_dist*max_epipolar_dist;

    const std::vector<MapLandmark *> & frame_1_3D_pts = frame_1->map_points_;
    const std::vector<MapLandmark *> & frame_2_3D_pts = frame_2->map_points_;

    bool bCam1_Calibrated = frame_1->getCamCalibrated();
    bool bCam2_Calibrated = frame_2->getCamCalibrated();

    const IntrinsicBase * cam_intrinsic_1 = frame_1->getCameraIntrinsics();
    const IntrinsicBase * cam_intrinsic_2 = frame_2->getCameraIntrinsics();


    double startTime = omp_get_wtime();

    matching::IndMatches putative_matches_1_2_idx;
    matchingWithCascadeHashing(frame_1,frame_2,putative_matches_1_2_idx, desc_ratio,max_desc_d);

    for (matching::IndMatches::iterator iter_pm = putative_matches_1_2_idx.begin(); iter_pm != putative_matches_1_2_idx.end();++iter_pm)
    {
      const IndexT feat_1_id = iter_pm->i_;
      const IndexT feat_2_id = iter_pm->j_;

      if (frame_1_3D_pts[feat_1_id] || frame_2_3D_pts[feat_2_id])
      {
        continue;
      }

      // Compute epipolar distance
      const Vec2 pt_1 = (bCam1_Calibrated ? frame_1->getFeaturePosition(feat_1_id) : cam_intrinsic_1->remove_disto(frame_1->getFeaturePosition(feat_1_id))).cast<double>();
      const Vec2 pt_2 = (bCam2_Calibrated ? frame_2->getFeaturePosition(feat_2_id) : cam_intrinsic_2->remove_disto(frame_2->getFeaturePosition(feat_2_id))).cast<double>();
      // Epipolar line on frame_2
      Vec3 x1_F = F_21 * pt_1.homogeneous() ; // Epipolar line on frame_2
      // Compute distance from epipolar line
      double dF_12 = x1_F.dot(pt_2.homogeneous());
      double d_pt_2_ep_line = (dF_12 * dF_12) /  (x1_F.head<2>().squaredNorm()); // square distance of pt_2 from epipolar line in frame_2

      if (d_pt_2_ep_line > max_epipolar_d_2)
      {
        continue;
      }

      map_matches_1_2_idx[feat_1_id] = feat_2_id;
    }

    std::cout<<"Matching: Epipolar 2D-2D ("<<omp_get_wtime() - startTime<<")\n";

  }






};

}
}
