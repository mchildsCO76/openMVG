
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <iterator>
#include <openMVG/vsslam/matching/Abstract_FeatureMatcher.hpp>
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include <openMVG/vsslam/tracking/PoseEstimation.hpp>

namespace openMVG  {
namespace VSSLAM  {

class Feat_Matcher_CascadeHashing : public Abstract_FeatureMatcher
{
private:

  void purgeCandidateMatches
  (
    Hash_Map<IndexT,IndexT> & map_putative_matches_1_2_idx,
    Hash_Map<IndexT,IndexT> & map_matches_1_2_idx
  )
  {
    // Check that two are not matching the same point
    bool bOkMatch = true;

    Hash_Map<IndexT,IndexT>::iterator iter_p_match_i = map_putative_matches_1_2_idx.begin();
    Hash_Map<IndexT,IndexT>::iterator iter_p_match_j;
    while (iter_p_match_i != map_putative_matches_1_2_idx.end())
    {
      iter_p_match_j = std::next(iter_p_match_i);
      while (iter_p_match_j != map_putative_matches_1_2_idx.end())
      {
        // Two features are mapped in the same feature
        if (iter_p_match_i->second == iter_p_match_j->second)
        {
          bOkMatch = false;
          iter_p_match_j = map_putative_matches_1_2_idx.end();
        }
        else
        {
          iter_p_match_j++;
        }
      }
      if (bOkMatch)
      {
        map_matches_1_2_idx[iter_p_match_i->first] = iter_p_match_i->second;
      }
      else
      {
        // reset duplicate flag
        bOkMatch = true;
      }
      iter_p_match_i++;
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
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  )
  {
    const std::string descriptor_type = frame_1->getRegions()->Type_id();

    if (descriptor_type == typeid(float).name())
    {
      matchingWithCascadeHashing_All_All_2D_2D<float>(frame_1,frame_2,vec_putative_matches,desc_ratio,feat_scale_ratio,max_desc_d);
    }
    else if (descriptor_type == typeid(double).name())
    {
      matchingWithCascadeHashing_All_All_2D_2D<double>(frame_1,frame_2,vec_putative_matches,desc_ratio,feat_scale_ratio,max_desc_d);
    }
    else if (descriptor_type == typeid(unsigned char).name())
    {
      matchingWithCascadeHashing_All_All_2D_2D<unsigned char>(frame_1,frame_2,vec_putative_matches,desc_ratio,feat_scale_ratio,max_desc_d);
    }
  }

  template <typename ScalarT>
  void matchingWithCascadeHashing_All_All_2D_2D
  (
    const Frame * frame_1,
    const Frame * frame_2,
    matching::IndMatches & vec_putative_matches,
    const float desc_ratio_sq = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d_sq = std::numeric_limits<float>::infinity()
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
      desc_ratio_sq);

    const std::vector<float> & frame_1_feat_scale = frame_1->pts_scale_;
    const std::vector<float> & frame_2_feat_scale = frame_2->pts_scale_;

    vec_putative_matches.reserve(vec_nn_ratio_idx.size());
    for (size_t k=0; k < vec_nn_ratio_idx.size(); ++k)
    {
      const size_t index = vec_nn_ratio_idx[k];
      // Check if distance of best descriptor is below threshold (useful?)
      if (pvec_distances[index*2] > max_desc_d_sq)
        continue;

      // Check if scale ratio between features is less than threshold
      const float f_scale_1 = frame_1_feat_scale[pvec_indices[index*2].j_];
      const float f_scale_2 = frame_2_feat_scale[pvec_indices[index*2].i_];
      if ((f_scale_1<f_scale_2 ? f_scale_2/f_scale_1 : f_scale_1/f_scale_2) > feat_scale_ratio)
        continue;

      vec_putative_matches.emplace_back(pvec_indices[index*2].j_, pvec_indices[index*2].i_);
    }
    // Remove duplicates
    matching::IndMatch::getDeduplicated(vec_putative_matches);

    // Remove matches that have the same (X,Y) coordinates
    matching::IndMatchDecorator<float> matchDeduplicator(vec_putative_matches,
        frame_1_regions->GetRegionsPositions(), frame_2_regions->GetRegionsPositions());
    matchDeduplicator.getDeduplicated(vec_putative_matches);


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
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float desc_ratio_2 = desc_ratio*desc_ratio;
    const float max_desc_d_2 = max_desc_d*max_desc_d;
    matchingWithCascadeHashing(frame_1,frame_2,vec_putative_matches_1_2_idx,desc_ratio_2,feat_scale_ratio,max_desc_d_2);
  }

  // Matching all available 3D points with the 2D features of frame_2
  void matching_AllAll_3D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const std::vector<MapLandmark *> & vec_3D_pts,
    const Frame * frame,
    Hash_Map<MapLandmark* ,IndexT> & matches_3D_pts_frame_idx,
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float desc_ratio_sq = desc_ratio*desc_ratio;
    const float max_desc_d_sq = max_desc_d*max_desc_d;

    // 3D data
    const size_t vec_3D_pts_size = vec_3D_pts.size();
    // Frame data
    features::Regions * const frame_regions = frame->getRegions();
    const std::vector<MapLandmark *> & frame_3D_pts = frame->map_points_;
    const std::vector<Vec2> & frame_pts = frame->pts_undist_;
    const IntrinsicBase * cam_intrinsic = frame->getCameraIntrinsics();

    // -------------------
    // -- Match features of prev_frame(that are already triangulated and in map) with the features of the new frame
    // -------------------
    std::vector<IndexT>matches_3D_frame_idx(vec_3D_pts_size,UndefinedIndexT);

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t pt_3D_i=0; pt_3D_i < vec_3D_pts_size; ++pt_3D_i)
    {
      MapLandmark * pt_3D = vec_3D_pts[pt_3D_i];

      // Check if its triangulated
      if (!pt_3D)
        continue;

      // Already matched
      if (matches_3D_pts_frame_idx.find(pt_3D) != matches_3D_pts_frame_idx.end())
        continue;

      // Project the point to frame 2
      Vec3 pt_3D_frame;
      PoseEstimator::getRelativePointPosition(pt_3D->X_,pt_3D->ref_frame_,pt_3D_frame,frame);
      if (pt_3D_frame(2) < 0)
        continue;

      // Project the point to image plane of frame 2
      const Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(cam_intrinsic->have_disto()?cam_intrinsic->add_disto(pt_3D_frame.hnormalized()):pt_3D_frame.hnormalized());

      // Check if projection is actually in the image borders
      if (!frame->isPointInFrame(pt_3D_frame_projected))
        continue;

      float pt_3D_scale = pt_3D->getFeatureMeanScale();
      
      // Raw descriptors
      void * pt_3D_raw_desc = pt_3D->feat_best_desc_;
      void * frame_pt_desc_raw;

      //TODO: Get possible candidates through grid

      // Loop through all the candidates and find the best and second best
      IndexT best_frame_feat_id = UndefinedIndexT;
      IndexT second_best_frame_feat_id = UndefinedIndexT;
      double best_desc_d_sq = max_desc_d_sq;
      double second_best_desc_d_sq = std::numeric_limits<double>::infinity();
      double desc_d_sq;

      for (size_t pt_f_i=0; pt_f_i < frame_pts.size(); ++pt_f_i)
      {
        // This point in the frame is already taken
        if (frame_3D_pts[pt_f_i])
          continue;

        // Check if points have feature scale ratio below threshold
        float frame_pt_scale = frame->getFeatureScale(pt_f_i);
        if ((frame_pt_scale<pt_3D_scale ? pt_3D_scale/frame_pt_scale : frame_pt_scale/pt_3D_scale) > feat_scale_ratio)
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_regions, pt_f_i, &frame_pt_desc_raw);

        // Compute distance_sq
        desc_d_sq = featureExtractor_->SquaredDescriptorDistance(pt_3D_raw_desc,frame_pt_desc_raw);

        // Save if in best two
        if (desc_d_sq < best_desc_d_sq)
        {
          // Former best move to second
          second_best_desc_d_sq = best_desc_d_sq;
          second_best_frame_feat_id = best_frame_feat_id;
          // Save as best
          best_desc_d_sq = desc_d_sq;
          best_frame_feat_id = pt_f_i;
        }
        else if (desc_d_sq < second_best_desc_d_sq)
        {
          // Save as second best
          second_best_desc_d_sq = desc_d_sq;
          second_best_frame_feat_id = pt_f_i;
        }
      }

      // Detect best match
      if (best_frame_feat_id != UndefinedIndexT)
      {
        if (second_best_frame_feat_id != UndefinedIndexT && (best_desc_d_sq / second_best_desc_d_sq) < desc_ratio_sq)
        {
          // Best is unique enough
          matches_3D_frame_idx[pt_3D_i] = best_frame_feat_id;
        }
        else
        {
          // Best is unique
          matches_3D_frame_idx[pt_3D_i] = best_frame_feat_id;
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
    const size_t max_px_d = 15, //20*20
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float desc_ratio_sq = desc_ratio*desc_ratio;
    const float max_desc_d_sq = max_desc_d * max_desc_d;
    const float max_px_d_sq = max_px_d * max_px_d;

    // 3D data
    const size_t vec_3D_pts_size = vec_3D_pts.size();
    // Frame data
    features::Regions * const frame_regions = frame->getRegions();
    const std::vector<MapLandmark *> & frame_3D_pts = frame->map_points_;
    const std::vector<Vec2> & frame_pts = frame->pts_undist_;
    const IntrinsicBase * cam_intrinsic = frame->getCameraIntrinsics();

    // -------------------
    // -- Match features of prev_frame(that are already triangulated and in map) with the features of the new frame
    // -------------------
    std::vector<IndexT>matches_3D_frame_idx(vec_3D_pts_size,UndefinedIndexT);

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t pt_3D_i=0; pt_3D_i < vec_3D_pts_size; ++pt_3D_i)
    {
      MapLandmark * pt_3D = vec_3D_pts[pt_3D_i];

      // Check if landmark exists
      if (!pt_3D)
        continue;

      // Already matched
      if (matches_3D_pts_frame_idx.find(pt_3D) != matches_3D_pts_frame_idx.end())
        continue;

      // Project the point to frame 2
      Vec3 pt_3D_frame;
      PoseEstimator::getRelativePointPosition(pt_3D->X_,pt_3D->ref_frame_,pt_3D_frame,frame);
      if (pt_3D_frame(2) < 0)
        continue;

      // Project the point to image plane of frame 2
      const Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(cam_intrinsic->have_disto()?cam_intrinsic->add_disto(pt_3D_frame.hnormalized()):pt_3D_frame.hnormalized());

      // Check if projection is actually in the image borders
      if (!frame->isPointInFrame(pt_3D_frame_projected))
        continue;

      float pt_3D_scale = pt_3D->getFeatureMeanScale();


      // Raw descriptors
      void * pt_3D_raw_desc = pt_3D->feat_best_desc_;
      void * frame_pt_desc_raw;

      //TODO: Get possible candidates through grid

      // Loop through all the candidates and find the best and second best
      IndexT best_frame_feat_id = UndefinedIndexT;
      IndexT second_best_frame_feat_id = UndefinedIndexT;
      double best_desc_d_sq = max_desc_d_sq;
      double second_best_desc_d_sq = std::numeric_limits<double>::infinity();
      double desc_d_sq;

      for (size_t pt_f_i=0; pt_f_i < frame_pts.size(); ++pt_f_i)
      {
        // This point in the frame is already taken
        if (frame_3D_pts[pt_f_i])
          continue;

        // Check if points are close enough
        if ((frame_pts[pt_f_i] - pt_3D_frame_projected).squaredNorm() > max_px_d_sq)
          continue;

        // Check if points have feature scale ratio below threshold
        float frame_pt_scale = frame->getFeatureScale(pt_f_i);
        if ((frame_pt_scale<pt_3D_scale ? pt_3D_scale/frame_pt_scale : frame_pt_scale/pt_3D_scale) > feat_scale_ratio)
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_regions, pt_f_i, &frame_pt_desc_raw);

        // Compute distance_sq
        desc_d_sq = featureExtractor_->SquaredDescriptorDistance(pt_3D_raw_desc,frame_pt_desc_raw);

        // Save if in best two
        if (desc_d_sq < best_desc_d_sq)
        {
          // Former best move to second
          second_best_desc_d_sq = best_desc_d_sq;
          second_best_frame_feat_id = best_frame_feat_id;
          // Save as best
          best_desc_d_sq = desc_d_sq;
          best_frame_feat_id = pt_f_i;
        }
        else if (desc_d_sq < second_best_desc_d_sq)
        {
          // Save as second best
          second_best_desc_d_sq = desc_d_sq;
          second_best_frame_feat_id = pt_f_i;
        }
      }

      // Detect best match
      if (best_frame_feat_id != UndefinedIndexT)
      {
        if (second_best_frame_feat_id != UndefinedIndexT && (best_desc_d_sq / second_best_desc_d_sq) < desc_ratio_sq)
        {
          // Best is unique enough
          matches_3D_frame_idx[pt_3D_i] = best_frame_feat_id;
        }
        else
        {
          // Best is unique
          matches_3D_frame_idx[pt_3D_i] = best_frame_feat_id;
        }
      }
    }

    // -------------------
    // -- Purge matches
    // -------------------
    purgeCandidateMatches(vec_3D_pts,matches_3D_frame_idx,matches_3D_pts_frame_idx);

  }


  // Match frame_2 to 3D points of frame_1
  void matching_Projection_3D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_idx,
    const size_t max_px_d = 15, //20*20
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float desc_ratio_sq = desc_ratio*desc_ratio;
    const float max_desc_d_sq = max_desc_d * max_desc_d;
    const float max_px_d_sq = max_px_d * max_px_d;

    // 3D data
    const std::vector<MapLandmark *> & frame_1_3D_pts = frame_1->map_points_;
    const std::vector<MapLandmark *> & frame_2_3D_pts = frame_1->map_points_;
    const size_t frame_1_3D_pts_size = frame_1_3D_pts.size();
    const size_t frame_2_3D_pts_size = frame_2_3D_pts.size();
    
    // Frame data
    features::Regions * const frame_2_regions = frame_2->getRegions();
    const std::vector<Vec2> & frame_2_2D_pts = frame_2->pts_undist_;
    const IntrinsicBase * cam_intrinsic_2 = frame_2->getCameraIntrinsics();

    // -------------------
    // -- Match features of prev_frame(that are already triangulated and in map) with the features of the new frame_2
    // -------------------
    std::vector<IndexT>matches_3D_frame_idx(frame_1_3D_pts_size,UndefinedIndexT);

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t pt_f1_i=0; pt_f1_i < frame_1_3D_pts_size; ++pt_f1_i)
    {
      MapLandmark * pt_3D = frame_1_3D_pts[pt_f1_i];

      // Check if landmark exists
      if (!pt_3D)
        continue;

      // Already matched
      if (matches_3D_pts_frame_idx.find(pt_3D) != matches_3D_pts_frame_idx.end())
        continue;

      // Project the point to frame_2 2
      Vec3 pt_3D_frame_2;
      PoseEstimator::getRelativePointPosition(pt_3D->X_,pt_3D->ref_frame_,pt_3D_frame_2,frame_2);
      if (pt_3D_frame_2(2) < 0)
        continue;

      // Project the point to image plane of frame_2 2
      const Vec2 pt_3D_frame_2_projected = cam_intrinsic_2->cam2ima(cam_intrinsic_2->have_disto()?cam_intrinsic_2->add_disto(pt_3D_frame_2.hnormalized()):pt_3D_frame_2.hnormalized());

      // Check if projection is actually in the image borders
      if (!frame_2->isPointInFrame(pt_3D_frame_2_projected))
        continue;

      float frame_1_pt_scale = frame_1->getFeatureScale(pt_f1_i);

      // Compute viewing angle between the point and normal which was last seen
      const Vec3 normal_X = (pt_3D->getWorldPosition() - frame_2->O_w_).normalized();
      const Vec3 & normal_last = pt_3D->getLastNormal();

      // Compute viewing angle
      float angle_factor = radiusByViewingAngle(normal_X.dot(normal_last));   // [1,5]
      // enlarge the area if the viewing angle is bigger
      float max_px_d_2 = (max_px_d * max_px_d) * angle_factor;

      // Raw descriptors
      void * pt_3D_raw_desc = pt_3D->feat_best_desc_;
      void * frame_2_pt_desc_raw;

      //TODO: Get possible candidates through grid

      // Loop through all the candidates and find the best and second best
      IndexT best_frame_feat_id = UndefinedIndexT;
      IndexT second_best_frame_feat_id = UndefinedIndexT;
      double best_desc_d_sq = max_desc_d_sq;
      double second_best_desc_d_sq = std::numeric_limits<double>::infinity();
      double desc_d_sq;

      for (size_t pt_f2_i=0; pt_f2_i < frame_2_3D_pts_size; ++pt_f2_i)
      {
        // This feature in the frame_2 is already taken
        if (frame_2->map_points_[pt_f2_i])
          continue;


        // Check if points are close enough
        if ((frame_2_2D_pts[pt_f2_i] - pt_3D_frame_2_projected).squaredNorm() > max_px_d_sq)
          continue;

        // Check if points have feature scale ratio below threshold
        float frame_2_pt_scale = frame_1->getFeatureScale(pt_f1_i);
        if ((frame_1_pt_scale<frame_2_pt_scale ? frame_2_pt_scale/frame_1_pt_scale : frame_1_pt_scale/frame_2_pt_scale) > feat_scale_ratio)
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_2_regions, pt_f2_i, &frame_2_pt_desc_raw);

        // Compute distance_2
        desc_d_sq = featureExtractor_->SquaredDescriptorDistance(pt_3D_raw_desc,frame_2_pt_desc_raw);

        // Save if in best two
        if (desc_d_sq < best_desc_d_sq)
        {
          // Former best move to second
          second_best_desc_d_sq = best_desc_d_sq;
          second_best_frame_feat_id = best_frame_feat_id;
          // Save as best
          best_desc_d_sq = desc_d_sq;
          best_frame_feat_id = pt_f2_i;
        }
        else if (desc_d_sq < second_best_desc_d_sq)
        {
          // Save as second best
          second_best_desc_d_sq = desc_d_sq;
          second_best_frame_feat_id = pt_f2_i;
        }
      }

      // Detect best match
      if (best_frame_feat_id != UndefinedIndexT)
      {
        if (second_best_frame_feat_id != UndefinedIndexT && (best_desc_d_sq / second_best_desc_d_sq) < desc_ratio_sq)
        {
          // Best is unique enough
          matches_3D_frame_idx[pt_f1_i] = best_frame_feat_id;
        }
        else
        {
          // Best is unique
          matches_3D_frame_idx[pt_f1_i] = best_frame_feat_id;
        }
      }
    }

    // -------------------
    // -- Purge matches
    // -------------------
    purgeCandidateMatches(frame_1_3D_pts,matches_3D_frame_idx,matches_3D_pts_frame_idx);

  }


  void matching_EpipolarLine_2D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    const Mat3 & F_21,
    Hash_Map<IndexT,IndexT> & map_matches_1_2_idx,
    const size_t epi_line_sigma = 4,
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float desc_ratio_sq = desc_ratio*desc_ratio;
    const float max_desc_d_sq = max_desc_d*max_desc_d;
    const float epi_line_thresh = (epi_line_sigma*epi_line_sigma) * 3.84;

    bool bCam1_Calibrated = frame_1->getCamCalibrated();
    bool bCam2_Calibrated = frame_2->getCamCalibrated();

    const IntrinsicBase * cam_intrinsic_1 = frame_1->getCameraIntrinsics();
    const IntrinsicBase * cam_intrinsic_2 = frame_2->getCameraIntrinsics();

    const std::vector<MapLandmark *> & frame_1_3D_pts = frame_1->map_points_;
    const std::vector<MapLandmark *> & frame_2_3D_pts = frame_2->map_points_;

    features::Regions * const frame_1_regions = frame_1->getRegions();
    features::Regions * const frame_2_regions = frame_2->getRegions();

    // Compute epipole in frame_2
    Vec3 frame_1_center = frame_1->getCameraCenter();

    // Project camera 1 center to frame 2
    Vec3 frame_1_center_frame_2;
    PoseEstimator::getRelativePointPosition(frame_1_center,nullptr,frame_1_center_frame_2,frame_2);

    // Project the point to image plane of frame 2
    const Vec2 epipole_frame_2 = cam_intrinsic_2->cam2ima(cam_intrinsic_2->have_disto()?cam_intrinsic_2->add_disto(frame_1_center_frame_2.hnormalized()):frame_1_center_frame_2.hnormalized());


    matching::IndMatches map_putative_matches_1_2_idx;
    matchingWithCascadeHashing(frame_1,frame_2,map_putative_matches_1_2_idx,desc_ratio_sq,feat_scale_ratio,max_desc_d_sq);

    for (matching::IndMatches::iterator iter_pm = map_putative_matches_1_2_idx.begin(); iter_pm != map_putative_matches_1_2_idx.end();++iter_pm)
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

      const float & frame_1_pt_scale = frame_1->getFeatureScale(feat_1_id);
      const float & frame_2_pt_scale = frame_2->getFeatureScale(feat_2_id);

      // Check if points have feature scale ratio below threshold
      if ((frame_1_pt_scale<frame_2_pt_scale ? frame_2_pt_scale/frame_1_pt_scale : frame_1_pt_scale/frame_2_pt_scale) > feat_scale_ratio)
        continue;


      // Check if point is not too close to epipole
      if ((epipole_frame_2 - pt_2).squaredNorm() < 100)
        continue;

      // Epipolar line on frame_2
      Vec3 x1_F = F_21 * pt_1.homogeneous() ; // Epipolar line on frame_2
      // Compute distance from epipolar line
      double dF_12 = x1_F.dot(pt_2.homogeneous());
      double d_pt_2_ep_line = (dF_12 * dF_12) /  (x1_F.head<2>().squaredNorm()); // square distance of pt_2 from epipolar line in frame_2

      if (d_pt_2_ep_line > epi_line_thresh)
      {
        continue;
      }

      map_matches_1_2_idx[feat_1_id] = feat_2_id;
    }


/*
    Hash_Map<IndexT,IndexT> map_putative_matches_1_2_idx;
    // For each point in frame 1 we try to find a match on epipolar line in step two
    for (size_t pt_f1_i = 0; pt_f1_i<frame_1_3D_pts.size(); ++pt_f1_i)
    {
      // If point is already matched we skip it
      if (frame_1_3D_pts[pt_f1_i])
        continue;

      const Vec2 & frame_1_pt = frame_1->getFeaturePosition(pt_f1_i);
      const float & frame_1_pt_scale = frame_1->getFeatureScale(pt_f1_i);

      // Raw descriptors
      void * frame_1_pt_desc_raw;
      void * frame_2_pt_desc_raw;

      // Get frame_1 descriptor
      featureExtractor_->getDescriptorRaw(frame_1_regions, pt_f1_i, &frame_1_pt_desc_raw);


      IndexT best_frame_feat_id = UndefinedIndexT;
      IndexT second_best_frame_feat_id = UndefinedIndexT;
      double best_desc_d_sq = max_desc_d_sq;
      double second_best_desc_d_sq = std::numeric_limits<double>::infinity();
      double desc_d_sq;

      for (size_t pt_f2_i=0; pt_f2_i < frame_2_3D_pts.size(); ++pt_f2_i)
      {
        // This point is already associated
        if (frame_2_3D_pts[pt_f2_i])
          continue;

        const Vec2 & frame_2_pt = frame_2->getFeaturePosition(pt_f2_i);
        const float & frame_2_pt_scale = frame_2->getFeatureScale(pt_f2_i);

        // Check if points have feature scale ratio below threshold
        if ((frame_1_pt_scale<frame_2_pt_scale ? frame_2_pt_scale/frame_1_pt_scale : frame_1_pt_scale/frame_2_pt_scale) > feat_scale_ratio)
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_2_regions, pt_f2_i, &frame_2_pt_desc_raw);

        // Compute distance_sq
        desc_d_sq = featureExtractor_->SquaredDescriptorDistance(frame_1_pt_desc_raw,frame_2_pt_desc_raw);

        if (desc_d_sq > second_best_desc_d_sq || desc_d_sq > max_desc_d_sq)
          continue;

        // Check if point is not too close to epipole
        if ((epipole_frame_2 - frame_2_pt).squaredNorm() < 100)
          continue;


        // Epipolar line on frame_2
        Vec3 F_x1 = F_21 * frame_1_pt.homogeneous();

        // Compute distance from epipolar line
        double err = F_x1.dot(frame_2_pt.homogeneous());
        double d_pt_2_ep_line_sq = (err * err) /  (F_x1.head<2>().squaredNorm()); // square distance of pt_2 from epipolar line in frame_2

        // Check if in validation area of epipolar line
        if (d_pt_2_ep_line_sq > epi_line_thresh)
          continue;

        // Save if in best two
        if (desc_d_sq < best_desc_d_sq)
        {
          // Former best move to second
          second_best_desc_d_sq = best_desc_d_sq;
          second_best_frame_feat_id = best_frame_feat_id;
          // Save as best
          best_desc_d_sq = desc_d_sq;
          best_frame_feat_id = pt_f2_i;
        }
        else if (desc_d_sq < second_best_desc_d_sq)
        {
          // Save as second best
          second_best_desc_d_sq = desc_d_sq;
          second_best_frame_feat_id = pt_f2_i;
        }
      }

      // Detect best match
      if (best_frame_feat_id != UndefinedIndexT)
      {
        if (second_best_frame_feat_id != UndefinedIndexT && (best_desc_d_sq / second_best_desc_d_sq) < desc_ratio_sq)
        {
          // Best is unique enough
          map_putative_matches_1_2_idx[pt_f1_i] = best_frame_feat_id;
        }
        else
        {
          // Best is unique
          map_putative_matches_1_2_idx[pt_f1_i] = best_frame_feat_id;
        }
      }

    }

    purgeCandidateMatches(map_putative_matches_1_2_idx,map_matches_1_2_idx);
*/

    //matchingWithCascadeHashing(frame_1,frame_2,vec_putative_matches_1_2_idx,desc_ratio_2,feat_scale_ratio,max_desc_d_2);
    //matchingWithCascadeHashing(frame_1,frame_2,putative_matches_1_2_idx, desc_ratio,std::numeric_limits<float>::infinity(),max_desc_d_sq);
    /*
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

      // Check if point is not too close to epipole
      if ((epipole_frame_2 - pt_2).squaredNorm() < 100)
        continue;

      // Epipolar line on frame_2
      Vec3 x1_F = F_21 * pt_1.homogeneous() ; // Epipolar line on frame_2
      // Compute distance from epipolar line
      double dF_12 = x1_F.dot(pt_2.homogeneous());
      double d_pt_2_ep_line = (dF_12 * dF_12) /  (x1_F.head<2>().squaredNorm()); // square distance of pt_2 from epipolar line in frame_2

      if (d_pt_2_ep_line > max_px_epi_d_2)
      {
        continue;
      }

      map_matches_1_2_idx[feat_1_id] = feat_2_id;
    }
     */
  }
};

}
}
