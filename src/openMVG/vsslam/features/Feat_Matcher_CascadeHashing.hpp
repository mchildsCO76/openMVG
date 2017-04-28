// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/types.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/matching_filters.hpp"
#include "openMVG/matching/matcher_cascade_hashing.hpp"

#include <openMVG/vsslam/system/Frame.hpp>
#include <openMVG/vsslam/vsslam_parameters.hpp>
#include <openMVG/vsslam/features/Abstract_Feature_Matcher.hpp>

namespace openMVG {
namespace vsslam {

class Feat_Matcher_CascadeHashing : public Abstract_Feature_Matcher
{
private:
  matching::CascadeHasher cascade_hasher_;

  void matchingWithCascadeHashing
    (
      const Frame * frame_1,
      const Frame * frame_2,
      matching::IndMatches & vec_putative_matches,
      const float f_desc_ratio_sq = 0.8,
      const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
      const float f_max_desc_d = std::numeric_limits<float>::infinity()
    )
    {
      const std::string descriptor_type = frame_1->getRegions()->Type_id();

      if (descriptor_type == typeid(float).name())
      {
        matchingWithCascadeHashing_All_All_2D_2D<float>(frame_1,frame_2,vec_putative_matches,f_desc_ratio_sq,f_feat_scale_ratio,f_max_desc_d);
      }
      else if (descriptor_type == typeid(double).name())
      {
        matchingWithCascadeHashing_All_All_2D_2D<double>(frame_1,frame_2,vec_putative_matches,f_desc_ratio_sq,f_feat_scale_ratio,f_max_desc_d);
      }
      else if (descriptor_type == typeid(unsigned char).name())
      {
        matchingWithCascadeHashing_All_All_2D_2D<unsigned char>(frame_1,frame_2,vec_putative_matches,f_desc_ratio_sq,f_feat_scale_ratio,f_max_desc_d);
      }
    }

  template <typename ScalarT>
  void matchingWithCascadeHashing_All_All_2D_2D
  (
    const Frame * frame_1,
    const Frame * frame_2,
    matching::IndMatches & vec_putative_matches,
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  )
  {
    double f_max_desc_d_sq = f_max_desc_d*f_max_desc_d;

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
      f_desc_ratio*f_desc_ratio);

    const std::vector<float> & frame_1_feat_scale = frame_1->getFeatureScaleVector();
    const std::vector<float> & frame_2_feat_scale = frame_2->getFeatureScaleVector();

    vec_putative_matches.reserve(vec_nn_ratio_idx.size());
    for (size_t k=0; k < vec_nn_ratio_idx.size(); ++k)
    {
      const size_t index = vec_nn_ratio_idx[k];
      // Check if distance of best descriptor is below threshold (useful?)
      if (pvec_distances[index*2] > f_max_desc_d_sq)
        continue;

      // Check if scale ratio between features is less than threshold
      const float f_scale_1 = frame_1_feat_scale[pvec_indices[index*2].j_];
      const float f_scale_2 = frame_2_feat_scale[pvec_indices[index*2].i_];
      if ((f_scale_1<f_scale_2 ? f_scale_2/f_scale_1 : f_scale_1/f_scale_2) > f_feat_scale_ratio)
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
  Feat_Matcher_CascadeHashing
  (
    std::shared_ptr<VSSLAM_Parameters> & params,
    Abstract_Feature_Extractor * feat_extractor
  ) : Abstract_Feature_Matcher(params)
  {
    cascade_hasher_.Init(feat_extractor->getDescriptorLength());
  }

  void purgeCandidateMatches
  (
    const std::vector<MapLandmark *> & map_landmarks,
    std::vector<IndexT> & putative_matches_3D_frame_idx,
    Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_idx
  )
  {
    // Check that two landmarks are not matching to the same landmark
    bool b_match_ok = true;

    for (size_t i=0; i<putative_matches_3D_frame_idx.size(); ++i)
    {
      if (putative_matches_3D_frame_idx[i] == UndefinedIndexT)
        continue;

      for (size_t j=i+1; j<putative_matches_3D_frame_idx.size(); ++j)
      {
        if (putative_matches_3D_frame_idx[j] == UndefinedIndexT)
          continue;
        // if value is repeated we delete both matches (matches have to be unique)
        if (putative_matches_3D_frame_idx[i] == putative_matches_3D_frame_idx[j])
        {
          b_match_ok = false;
          putative_matches_3D_frame_idx[j] = UndefinedIndexT;
        }
      }

      if (b_match_ok)
      {
        matches_3D_pts_frame_idx[map_landmarks[i]] = putative_matches_3D_frame_idx[i];
      }
      else
      {
        // reset flag
        b_match_ok = true;
      }
    }


  }


  void matching_AllAll_2D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor,
    const Frame * frame_1,
    const Frame * frame_2,
    matching::IndMatches & vec_putative_matches_1_2_idx,
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float f_desc_ratio_sq = f_desc_ratio*f_desc_ratio;
    const float f_max_desc_d_sq = f_max_desc_d*f_max_desc_d;
    matchingWithCascadeHashing(frame_1,frame_2,vec_putative_matches_1_2_idx,f_desc_ratio_sq,f_feat_scale_ratio,f_max_desc_d_sq);
  }

  void matching_AllAll_3D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor,
    const std::vector<MapLandmark *> & vec_landmarks,
    const Frame * frame,
    Hash_Map<MapLandmark* ,IndexT> & matches_landmark_frame_idx,
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float f_desc_ratio_sq = f_desc_ratio*f_desc_ratio;
    const float f_max_desc_d_sq = f_max_desc_d * f_max_desc_d;

    // 3D data
    const size_t n_landmarks = vec_landmarks.size();
    const size_t & n_feats_frame = frame->getNumberOfFeatures();

    // Frame data
    features::Regions * const frame_regions = frame->getRegions();
    const std::vector<MapLandmark *> & vec_landmarks_frame = frame->getLandmarks();
    const std::vector<Vec2> & frame_2D_pts = frame->getFeaturePointVector();

    // -------------------
    // -- Match landmarks with features in frame
    // -------------------
    std::vector<IndexT> vec_putative_matches_landmark_frame_idx(n_feats_frame,UndefinedIndexT);

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (IndexT pt_landmark_i=0; pt_landmark_i < n_landmarks; ++pt_landmark_i)
    {
      MapLandmark * map_landmark = vec_landmarks[pt_landmark_i];


      // Check if the point is already matched
      if (matches_landmark_frame_idx.find(map_landmark) != matches_landmark_frame_idx.end())
        continue;

      // Project to current frame
      Vec2 pt_2D_frame_2_projected;
      if (!frame->getProjectedPoint(map_landmark,pt_2D_frame_2_projected))
        continue;

      float map_landmark_scale = map_landmark->getFeatureMeanScale();

      // Raw descriptors
      void * map_landmark_raw_desc = map_landmark->getBestDesc();
      void * pt_2D_desc_frame;

      //TODO: Get possible candidates through grid

      // ------------------------------
      // -- Loop through all the candidates and find the best and second best
      // ------------------------------
      IndexT feat_id_best = UndefinedIndexT;
      IndexT feat_id_second = UndefinedIndexT;
      double f_desc_d_best_sq = f_max_desc_d_sq;
      double f_desc_d_second_sq = std::numeric_limits<double>::infinity();
      double f_desc_d_sq;

      // Loop through fetures in current frame
      for (IndexT pt_frame_i=0; pt_frame_i < n_feats_frame; ++pt_frame_i)
      {
        if (vec_landmarks_frame[pt_frame_i])
          continue;

        // Check if points have feature scale ratio below threshold
        float pt_2d_frame_scale = frame->getFeatureScale(pt_frame_i);
        if ((pt_2d_frame_scale<map_landmark_scale ? map_landmark_scale/pt_2d_frame_scale : pt_2d_frame_scale/map_landmark_scale) > f_feat_scale_ratio)
          continue;

        // Get candidate descriptor
        feature_extractor->getDescriptorRaw(frame_regions, pt_frame_i, &pt_2D_desc_frame);

        // Compute distance_sq
        f_desc_d_sq = feature_extractor->getDescriptorDistanceSquared(map_landmark_raw_desc,pt_2D_desc_frame);

        // Save if in best two
        if (f_desc_d_sq < f_desc_d_best_sq)
        {
          // Former best move to second
          feat_id_second = feat_id_best;
          f_desc_d_second_sq = f_desc_d_best_sq;
          // Save as best
          feat_id_best = pt_frame_i;
          f_desc_d_best_sq = f_desc_d_sq;
        }
        else if (f_desc_d_sq < f_desc_d_second_sq)
        {
          // Save as second best
          feat_id_second = pt_frame_i;
          f_desc_d_second_sq = f_desc_d_sq;
        }
      }

      // Detect best match
      if (feat_id_best != UndefinedIndexT)
      {
        // If we have to matches we check ratio otherwise directly select best
        if (feat_id_second != UndefinedIndexT)
        {
          if ((f_desc_d_best_sq / f_desc_d_second_sq) < f_desc_ratio_sq)
            vec_putative_matches_landmark_frame_idx[pt_landmark_i] = feat_id_best;
        }
        else
        {
          vec_putative_matches_landmark_frame_idx[pt_landmark_i] = feat_id_best;
        }
      }
    }
    // -------------------
    // -- Purge matches - check if two map_landmarks match the same feature
    // -------------------
    purgeCandidateMatches(vec_landmarks_frame,vec_putative_matches_landmark_frame_idx,matches_landmark_frame_idx);

  }

  void matching_Projection_3D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor,
    const Frame * frame_1,
    const Frame * frame_2,
    Hash_Map<MapLandmark *,IndexT> & matches_landmark_frame_idx,
    const float f_max_px_d = 15,
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float f_desc_ratio_sq = f_desc_ratio*f_desc_ratio;
    const float f_max_desc_d_sq = f_max_desc_d * f_max_desc_d;
    const float f_max_px_d_sq = f_max_px_d * f_max_px_d;

    // 3D data
    const std::vector<MapLandmark *> & vec_landmarks_frame_1 = frame_1->getLandmarks();
    const std::vector<MapLandmark *> & vec_landmarks_frame_2 = frame_2->getLandmarks();
    const size_t & n_feats_frame_1 = frame_1->getNumberOfFeatures();
    const size_t & n_feats_frame_2 = frame_2->getNumberOfFeatures();

    // Frame data
    features::Regions * const frame_2_regions = frame_2->getRegions();
    const std::vector<Vec2> & frame_2_2D_pts = frame_2->getFeaturePointVector();

    // -------------------
    // -- Try to match features tracked in prev frame by projection to current and looking in the neighborhood
    // -------------------
    std::vector<IndexT> vec_putative_matches_frame_1_2_idx(n_feats_frame_1,UndefinedIndexT);

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (IndexT pt_frame_1_i=0; pt_frame_1_i < n_feats_frame_1; ++pt_frame_1_i)
    {
      MapLandmark * map_landmark = vec_landmarks_frame_1[pt_frame_1_i];

      if (!map_landmark)
        continue;

      // Check if the point is already matched
      if (matches_landmark_frame_idx.find(map_landmark) != matches_landmark_frame_idx.end())
        continue;

      // Project to current frame
      Vec2 pt_2D_frame_2_projected;
      if (!frame_2->getProjectedPoint(map_landmark,pt_2D_frame_2_projected))
        continue;

      float map_landmark_scale = map_landmark->getFeatureMeanScale();

      // Raw descriptors
      void * map_landmark_raw_desc = map_landmark->getBestDesc();
      void * pt_2D_desc_frame;

      //TODO: Get possible candidates through grid

      // ------------------------------
      // -- Loop through all the candidates and find the best and second best
      // ------------------------------
      IndexT feat_id_best = UndefinedIndexT;
      IndexT feat_id_second = UndefinedIndexT;
      double f_desc_d_best_sq = f_max_desc_d_sq;
      double f_desc_d_second_sq = std::numeric_limits<double>::infinity();
      double f_desc_d_sq;

      // Loop through fetures in current frame
      for (IndexT pt_frame_2_i=0; pt_frame_2_i < n_feats_frame_2; ++pt_frame_2_i)
      {
        if (vec_landmarks_frame_2[pt_frame_2_i])
          continue;

        // Check if points are close enough
        if ((frame_2_2D_pts[pt_frame_2_i] - pt_2D_frame_2_projected).squaredNorm() > f_max_px_d_sq)
          continue;

        // Check if points have feature scale ratio below threshold
        float pt_2d_frame_scale = frame_2->getFeatureScale(pt_frame_2_i);
        if ((pt_2d_frame_scale<map_landmark_scale ? map_landmark_scale/pt_2d_frame_scale : pt_2d_frame_scale/map_landmark_scale) > f_feat_scale_ratio)
          continue;

        // Get candidate descriptor
        feature_extractor->getDescriptorRaw(frame_2_regions, pt_frame_2_i, &pt_2D_desc_frame);

        // Compute distance_sq
        f_desc_d_sq = feature_extractor->getDescriptorDistanceSquared(map_landmark_raw_desc,pt_2D_desc_frame);

        // Save if in best two
        if (f_desc_d_sq < f_desc_d_best_sq)
        {
          // Former best move to second
          feat_id_second = feat_id_best;
          f_desc_d_second_sq = f_desc_d_best_sq;
          // Save as best
          feat_id_best = pt_frame_2_i;
          f_desc_d_best_sq = f_desc_d_sq;
        }
        else if (f_desc_d_sq < f_desc_d_second_sq)
        {
          // Save as second best
          feat_id_second = pt_frame_2_i;
          f_desc_d_second_sq = f_desc_d_sq;
        }
      }

      // Detect best match
      if (feat_id_best != UndefinedIndexT)
      {
        // If we have to matches we check ratio otherwise directly select best
        if (feat_id_second != UndefinedIndexT)
        {
          if ((f_desc_d_best_sq / f_desc_d_second_sq) < f_desc_ratio_sq)
            vec_putative_matches_frame_1_2_idx[pt_frame_1_i] = feat_id_best;
        }
        else
        {
          vec_putative_matches_frame_1_2_idx[pt_frame_1_i] = feat_id_best;
        }
      }
    }

    // -------------------
    // -- Purge matches - check if two map_landmarks match the same feature
    // -------------------
    purgeCandidateMatches(vec_landmarks_frame_1,vec_putative_matches_frame_1_2_idx,matches_landmark_frame_idx);
  }

  void matching_Projection_3D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor,
    const std::vector<MapLandmark *> & vec_landmarks,
    const Frame * frame,
    Hash_Map<MapLandmark* ,IndexT> & matches_landmark_frame_idx,
    const float f_max_px_d = 15,
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float f_desc_ratio_sq = f_desc_ratio*f_desc_ratio;
    const float f_max_desc_d_sq = f_max_desc_d * f_max_desc_d;
    const float f_max_px_d_sq = f_max_px_d * f_max_px_d;

    // 3D data
    const size_t n_landmarks = vec_landmarks.size();
    const size_t & n_feats_frame = frame->getNumberOfFeatures();

    // Frame data
    features::Regions * const frame_regions = frame->getRegions();
    const std::vector<MapLandmark *> & vec_landmarks_frame = frame->getLandmarks();
    const std::vector<Vec2> & frame_2D_pts = frame->getFeaturePointVector();

    // -------------------
    // -- Match landmarks with features in frame
    // -------------------
    std::vector<IndexT> vec_putative_matches_landmark_frame_idx(n_landmarks,UndefinedIndexT);

    //#ifdef OPENMVG_USE_OPENMP
    //#pragma omp parallel for schedule(dynamic)
    //#endif
    for (IndexT pt_landmark_i=0; pt_landmark_i < n_landmarks; ++pt_landmark_i)
    {
      MapLandmark * map_landmark = vec_landmarks[pt_landmark_i];

      if (!map_landmark)
      {
        std::cout<<"LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL\n";
      }

      // Check if the point is already matched
      if (matches_landmark_frame_idx.find(map_landmark) != matches_landmark_frame_idx.end())
        continue;

      // Project to current frame
      Vec2 pt_2D_frame_projected;
      if (!frame->getProjectedPoint(map_landmark,pt_2D_frame_projected))
        continue;

      float map_landmark_scale = map_landmark->getFeatureMeanScale();

      // Raw descriptors
      void * map_landmark_raw_desc = map_landmark->getBestDesc();
      void * pt_2D_desc_frame;

      //TODO: Get possible candidates through grid

      // ------------------------------
      // -- Loop through all the candidates and find the best and second best
      // ------------------------------
      IndexT feat_id_best = UndefinedIndexT;
      IndexT feat_id_second = UndefinedIndexT;
      double f_desc_d_best_sq = f_max_desc_d_sq;
      double f_desc_d_second_sq = std::numeric_limits<double>::infinity();
      double f_desc_d_sq;

      // Loop through fetures in current frame
      for (IndexT pt_frame_i=0; pt_frame_i < n_feats_frame; ++pt_frame_i)
      {
        if (vec_landmarks_frame[pt_frame_i])
          continue;

        // Check if points are close enough
        if ((frame_2D_pts[pt_frame_i] - pt_2D_frame_projected).squaredNorm() > f_max_px_d_sq)
          continue;


        // Check if points have feature scale ratio below threshold
        float pt_2d_frame_scale = frame->getFeatureScale(pt_frame_i);
        if ((pt_2d_frame_scale<map_landmark_scale ? map_landmark_scale/pt_2d_frame_scale : pt_2d_frame_scale/map_landmark_scale) > f_feat_scale_ratio)
          continue;


        // Get candidate descriptor
        feature_extractor->getDescriptorRaw(frame_regions, pt_frame_i, &pt_2D_desc_frame);

        // Compute distance_sq
        f_desc_d_sq = feature_extractor->getDescriptorDistanceSquared(map_landmark_raw_desc,pt_2D_desc_frame);


        // Save if in best two
        if (f_desc_d_sq < f_desc_d_best_sq)
        {

          // Former best move to second
          feat_id_second = feat_id_best;
          f_desc_d_second_sq = f_desc_d_best_sq;
          // Save as best
          feat_id_best = pt_frame_i;
          f_desc_d_best_sq = f_desc_d_sq;
        }
        else if (f_desc_d_sq < f_desc_d_second_sq)
        {
          // Save as second best
          feat_id_second = pt_frame_i;
          f_desc_d_second_sq = f_desc_d_sq;
        }
      }

      // Detect best match
      if (feat_id_best != UndefinedIndexT)
      {
        // If we have to matches we check ratio otherwise directly select best
        if (feat_id_second != UndefinedIndexT)
        {
          if ((f_desc_d_best_sq / f_desc_d_second_sq) < f_desc_ratio_sq)
            vec_putative_matches_landmark_frame_idx[pt_landmark_i] = feat_id_best;
        }
        else
        {
          vec_putative_matches_landmark_frame_idx[pt_landmark_i] = feat_id_best;
        }
      }
    }

    // -------------------
    // -- Purge matches - check if two map_landmarks match the same feature
    // -------------------
    purgeCandidateMatches(vec_landmarks,vec_putative_matches_landmark_frame_idx,matches_landmark_frame_idx);

  }


  void matching_Epipolar_2D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    const Mat3 & F_21,
    Hash_Map<IndexT,IndexT> & map_matches_1_2_idx,
    const float f_max_px_d = 4,
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float f_desc_ratio_sq = f_desc_ratio*f_desc_ratio;
    const float f_max_desc_d_sq = f_max_desc_d * f_max_desc_d;
    const float f_max_px_d_sq = f_max_px_d * f_max_px_d;

    // 3D data
    const std::vector<MapLandmark *> & vec_landmarks_frame_1 = frame_1->getLandmarks();
    const std::vector<MapLandmark *> & vec_landmarks_frame_2 = frame_2->getLandmarks();
    const size_t & n_feats_frame_1 = frame_1->getNumberOfFeatures();
    const size_t & n_feats_frame_2 = frame_2->getNumberOfFeatures();

    // Frame data
    features::Regions * const frame_1_regions = frame_1->getRegions();
    const std::vector<Vec2> & frame_1_2D_pts = frame_1->getFeaturePointVector();
    features::Regions * const frame_2_regions = frame_2->getRegions();
    const std::vector<Vec2> & frame_2_2D_pts = frame_2->getFeaturePointVector();

    // Compute epipole in frame_2
    Vec2 epipole_1_frame_2;
    frame_2->getProjectedPoint(frame_1->getCameraCenter(), nullptr, epipole_1_frame_2);

    matching::IndMatches map_putative_matches_1_2_idx;
    matchingWithCascadeHashing(frame_1,frame_2,map_putative_matches_1_2_idx,f_desc_ratio_sq,f_feat_scale_ratio,f_max_desc_d_sq);


    for (matching::IndMatches::iterator match_it = map_putative_matches_1_2_idx.begin(); match_it != map_putative_matches_1_2_idx.end();++match_it)
    {
      const IndexT feat_id_1 = match_it->i_;
      const IndexT feat_id_2 = match_it->j_;

      // If either of features is already assigned we skip the match
      if (vec_landmarks_frame_1[feat_id_1] || vec_landmarks_frame_2[feat_id_2])
        continue;

      // Check feature scale compatibility
      float pt_2D_frame_1_scale = frame_1->getFeatureScale(feat_id_1);
      float pt_2D_frame_2_scale = frame_2->getFeatureScale(feat_id_2);

      if (pt_2D_frame_1_scale > 4 || pt_2D_frame_2_scale > 4)
        continue;

      if ((pt_2D_frame_1_scale<pt_2D_frame_2_scale ? pt_2D_frame_2_scale/pt_2D_frame_1_scale : pt_2D_frame_1_scale/pt_2D_frame_2_scale) > f_feat_scale_ratio)
        continue;

      // Get feature position
      Vec2 pt_2D_frame_1 = frame_1->getFeaturePosition(feat_id_1);
      Vec2 pt_2D_frame_2 = frame_2->getFeaturePosition(feat_id_2);

      // Check if its not too close to epipole
      if ((epipole_1_frame_2 - pt_2D_frame_2).squaredNorm() < 900 )
        continue;

      // Compute distance to epipolar line on frame 2
      // Epipolar line on frame_2
      Vec3 F_x1 = F_21 * pt_2D_frame_1.homogeneous() ;
      double x2_F_x1 = F_x1.dot(pt_2D_frame_2.homogeneous());
      // Compute distance from epipolar line
      double pt2_F_d = (x2_F_x1 * x2_F_x1) /  (F_x1.head<2>().squaredNorm()); // square distance of pt_2 from epipolar line in frame_2

      if (pt2_F_d > f_max_px_d_sq)
        continue;

      map_matches_1_2_idx[feat_id_1] = feat_id_2;
    }


  }


  /*
  void matching_Projection_3D_2D
  (
    const Abstract_Feature_Extractor * feature_extractor_,
    const Frame * frame,
    const std::vector<MapLandmark *> & vec_3D_pts,
    Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_idx,
    const float f_max_px_d = 15, //20*20
    const float f_desc_ratio = 0.8,
    const float f_feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float f_max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    const float f_desc_ratio_sq = f_desc_ratio*f_desc_ratio;
    const float f_max_desc_d_sq = f_max_desc_d * f_max_desc_d;
    const float f_max_px_d_sq = f_max_px_d * f_max_px_d;

    // 3D data
    const size_t vec_3D_pts_size = vec_3D_pts.size();
    // Frame data
    features::Regions * const frame_regions = frame->getRegions();
    const std::vector<MapLandmark *> & frame_3D_pts = frame->getLandmarks();
    const std::vector<Vec2> & frame_pts_2D = frame->getFeaturePointVector();
    const IntrinsicBase * cam_intrinsic = frame->getCameraIntrinsics();

    // -------------------
    // -- Try to match features tracked in prev frame by projection to current and looking in the neighborhood
    // -------------------
    std::vector<IndexT>matches_3D_frame_idx(vec_3D_pts_size,UndefinedIndexT);

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t m_l_i=0; m_l_i < vec_3D_pts_size; ++m_l_i)
    {
      MapLandmark * map_landmark = vec_3D_pts[m_l_i];

      if (!map_landmark)
        continue;

      // Check if the point is already matched
      if (matches_3D_pts_frame_idx.find(map_landmark) != matches_3D_pts_frame_idx.end())
        continue;

      // Project to current frame
      Vec2 pt_3D_frame_projected;
      if (!frame->getProjectedPointMapLandmark(map_landmark,pt_3D_frame_projected))
        continue;

      float map_landmark_scale = map_landmark->getFeatureMeanScale();

      // Raw descriptors
      void * map_landmark_raw_desc = map_landmark->getBestDesc();
      void * pt_2D_desc_frame;

      //TODO: Get possible candidates through grid

      // ------------------------------
      // -- Loop through all the candidates and find the best and second best
      // ------------------------------
      IndexT feat_id_best = UndefinedIndexT;
      IndexT feat_id_second = UndefinedIndexT;
      double f_desc_d_best_sq = f_max_desc_d_sq;
      double f_desc_d_second_sq = std::numeric_limits<double>::infinity();
      double f_desc_d_sq;

      // Loop through fetures in current frame
      for (size_t pt_i_frame=0; pt_i_frame < frame_pts_2D.size(); ++pt_i_frame)
      {
        if (frame_3D_pts[pt_i_frame])
          continue;

        // Check if points are close enough
        if ((frame_pts_2D[pt_i_frame] - pt_3D_frame_projected).squaredNorm() > f_max_px_d_sq)
          continue;

        // Check if points have feature scale ratio below threshold
        float pt_2d_frame_scale = frame->getFeatureScale(pt_i_frame);
        if ((pt_2d_frame_scale<map_landmark_scale ? map_landmark_scale/pt_2d_frame_scale : pt_2d_frame_scale/map_landmark_scale) > f_feat_scale_ratio)
          continue;

        // Get candidate descriptor
        feature_extractor_->getDescriptorRaw(frame_regions, pt_i_frame, &pt_2D_desc_frame);

        // Compute distance_sq
        f_desc_d_sq = feature_extractor_->SquaredDescriptorDistance(map_landmark_raw_desc,pt_2D_desc_frame);

        // Save if in best two
        if (f_desc_d_sq < f_desc_d_best_sq)
        {
          // Former best move to second
          feat_id_second = feat_id_best;
          f_desc_d_second_sq = f_desc_d_best_sq;
          // Save as best
          feat_id_best = pt_i_frame;
          f_desc_d_best_sq = f_desc_d_sq;
        }
        else if (f_desc_d_sq < f_desc_d_second_sq)
        {
          // Save as second best
          feat_id_second = pt_i_frame;
          f_desc_d_second_sq = f_desc_d_sq;
        }
      }

      // Detect best match
      if (feat_id_best != UndefinedIndexT)
      {
        // If we have to matches we check ratio otherwise directly select best
        if (feat_id_second != UndefinedIndexT)
        {
          if ((f_desc_d_best_sq / f_desc_d_second_sq) < f_desc_ratio_sq)
            matches_3D_frame_idx[m_l_i] = feat_id_best;
        }
        else
        {
          matches_3D_frame_idx[m_l_i] = feat_id_best;
        }
      }
    }

    // -------------------
    // -- Purge matches - check if two map_landmarks match the same feature
    // -------------------
    purgeCandidateMatches(frame_3D_pts,matches_3D_frame_idx,matches_3D_pts_frame_idx);
  }
  */


};

}
}
