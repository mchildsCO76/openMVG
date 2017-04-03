
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/vsslam/matching/Abstract_FeatureMatcher.hpp>
#include <openMVG/vsslam/tracking/PoseEstimation.hpp>

namespace openMVG  {
namespace VSSLAM  {

class Feat_Matcher_Basic : public Abstract_FeatureMatcher
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
    std::vector<IndexT> & vec_putative_matches_1_2_idx,
    Hash_Map<IndexT,IndexT> & map_matches_1_2_idx
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
          map_matches_1_2_idx[i] = vec_putative_matches_1_2_idx[i];
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

public:

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

    // Frame_1 data
    features::Regions * const frame_1_regions = frame_1->getRegions();
    const size_t frame_1_pts_size = frame_1->getNumberOfFeatures();
    // Frame_2 data
    features::Regions * const frame_2_regions = frame_2->getRegions();
    const size_t frame_2_pts_size = frame_2->getNumberOfFeatures();

    // -------------------
    // -- Match reference-current candidates
    // -------------------
    std::vector<IndexT>putative_matches_1_2_idx(frame_1_pts_size,UndefinedIndexT);


    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (IndexT p_i=0; p_i < frame_1_pts_size; ++p_i)
    {
      IndexT best_idx = UndefinedIndexT;
      IndexT second_best_idx = UndefinedIndexT;
      double best_distance_2 = max_desc_d_2;
      double second_best_distance_2 = std::numeric_limits<double>::infinity();
      double distance_2;

      // Raw descriptors
      void * ref_pt_desc_raw;
      void * candidate_pt_desc_raw;

      // Get descriptor of reference point
      featureExtractor_->getDescriptorRaw(frame_1_regions, p_i, &ref_pt_desc_raw);

      for (IndexT c_i=0; c_i < frame_2_pts_size; ++c_i)
      {
        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_2_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance
        distance_2 = featureExtractor_->SquaredDescriptorDistance(ref_pt_desc_raw,candidate_pt_desc_raw);

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

      const std::vector<float> & frame_1_feat_scale = frame_1->pts_scale_;
      const std::vector<float> & frame_2_feat_scale = frame_2->pts_scale_;

      // Detect best match
      if (best_idx != UndefinedIndexT && best_distance_2 < max_desc_d_2)
      {
        if (second_best_idx != UndefinedIndexT)
        {
          if ((best_distance_2 / second_best_distance_2) < desc_ratio_2)
          {

            // Check if scale ratio between features is less than threshold
            const float f_scale_1 = frame_1_feat_scale[p_i];
            const float f_scale_2 = frame_2_feat_scale[best_idx];
            if ((f_scale_1<f_scale_2 ? f_scale_2/f_scale_1 : f_scale_1/f_scale_2) < feat_scale_ratio)
              // Best is unique enough
              putative_matches_1_2_idx[p_i] = best_idx;
          }
        }
        else
        {
          const float f_scale_1 = frame_1_feat_scale[p_i];
          const float f_scale_2 = frame_2_feat_scale[best_idx];
          if ((f_scale_1<f_scale_2 ? f_scale_2/f_scale_1 : f_scale_1/f_scale_2) < feat_scale_ratio)
            // Best is unique enough
            putative_matches_1_2_idx[p_i] = best_idx;

        }

      }
    }
    // -------------------
    // -- Purge matches
    // -------------------
    purgeCandidateMatches(putative_matches_1_2_idx, vec_putative_matches_1_2_idx);
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
    const float max_desc_d_2 = std::numeric_limits<float>::infinity()
  ) override
  {
    const float desc_ratio_2 = desc_ratio*desc_ratio;

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

    double startTime = omp_get_wtime();

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
      PoseEstimator::getRelativePointPosition(pt_3D->X_,pt_3D->ref_frame_,pt_3D_frame,frame);
      if (pt_3D_frame(2) < 0)
        continue;

      // Project the point to image plane of frame 2
      const Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(cam_intrinsic->have_disto()?cam_intrinsic->add_disto(pt_3D_frame.hnormalized()):pt_3D_frame.hnormalized());

      // Check if projection is actually in the image borders
      if (!frame->isPointInFrame(pt_3D_frame_projected))
        continue;

      // Raw descriptors
      void * pt_3D_raw_desc = pt_3D->feat_best_desc_;
      void * candidate_pt_desc_raw;

      //TODO: Get possible candidates through grid
      IndexT best_idx = UndefinedIndexT;
      IndexT second_best_idx = UndefinedIndexT;
      double best_distance_2 = std::numeric_limits<double>::infinity();
      double second_best_distance_2 = std::numeric_limits<double>::infinity();
      double distance_2;

      for (size_t c_i=0; c_i < frame_pts.size(); ++c_i)
      {
        // This point in the frame is already taken
        if (frame->map_points_[c_i])
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance
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
      if (best_idx != UndefinedIndexT && best_distance_2 < max_desc_d_2)
      {
        if (second_best_idx != UndefinedIndexT  && second_best_distance_2 < max_desc_d_2)
        {
          if ((best_distance_2 / second_best_distance_2) < desc_ratio_2)
          {
            // Best is unique enough
            matches_3D_frame_idx[p_i] = best_idx;
          }
        }
        else
        {
          // Best is unique
          matches_3D_frame_idx[p_i] = best_idx;
        }
      }
    }

    std::cout<<"Matching: Projection 3D-2D ("<<omp_get_wtime() - startTime<<")\n";

    // -------------------
    // -- Purge matches
    // -------------------
    startTime = omp_get_wtime();

    purgeCandidateMatches(vec_3D_pts,matches_3D_frame_idx,matches_3D_pts_frame_idx);

    std::cout<<"Matching: Purge ("<<omp_get_wtime() - startTime<<")\n";

  }

  void matching_Projection_3D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame,
    const std::vector<MapLandmark *> & vec_3D_pts,
    Hash_Map<MapLandmark *,IndexT> & matches_3D_pts_frame_idx,
    const size_t win_size = 400, //20*20
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d_2 = std::numeric_limits<float>::infinity()
  ) override
  {
    const float desc_ratio_2 = desc_ratio*desc_ratio;
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

    double startTime = omp_get_wtime();

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
      PoseEstimator::getRelativePointPosition(pt_3D->X_,pt_3D->ref_frame_,pt_3D_frame,frame);
      if (pt_3D_frame(2) < 0)
        continue;

      // Project the point to image plane of frame 2
      const Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(cam_intrinsic->have_disto()?cam_intrinsic->add_disto(pt_3D_frame.hnormalized()):pt_3D_frame.hnormalized());

      // Check if projection is actually in the image borders
      if (!frame->isPointInFrame(pt_3D_frame_projected))
        continue;

      // Raw descriptors
      void * pt_3D_raw_desc = pt_3D->feat_best_desc_;
      void * candidate_pt_desc_raw;

      //TODO: Get possible candidates through grid

      IndexT best_idx = UndefinedIndexT;
      IndexT second_best_idx = UndefinedIndexT;
      double best_distance_2 = std::numeric_limits<double>::infinity();
      double second_best_distance_2 = std::numeric_limits<double>::infinity();
      double distance_2;

      for (size_t c_i=0; c_i < frame_pts.size(); ++c_i)
      {
        // This point in the frame is already taken
        if (frame->map_points_[c_i])
          continue;

        // Check if points are close enough
        if ((frame_pts[c_i] - pt_3D_frame_projected).squaredNorm() > win_size)
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance
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
      if (best_idx != UndefinedIndexT && best_distance_2 < max_desc_d_2)
      {
        if (second_best_idx != UndefinedIndexT  && second_best_distance_2 < max_desc_d_2)
        {
          if ((best_distance_2 / second_best_distance_2) < desc_ratio_2)
          {
            // Best is unique enough
            matches_3D_frame_idx[p_i] = best_idx;
          }
        }
        else
        {
          // Best is unique
          matches_3D_frame_idx[p_i] = best_idx;
        }
      }
    }

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Matching: Projection 3D-2D ("<<secsElapsed<<")\n";

    // -------------------
    // -- Purge matches
    // -------------------
    startTime = omp_get_wtime();

    purgeCandidateMatches(vec_3D_pts,matches_3D_frame_idx,matches_3D_pts_frame_idx);

    std::cout<<"Matching: Purge ("<<omp_get_wtime() - startTime<<")\n";
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
  }

  void matching_EpipolarLine_2D_2D
  (
    const Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    const Mat3 & F_21,
    Hash_Map<IndexT,IndexT> & map_matches_1_2_idx,
    const size_t max_px_epi_d = 4,
    const float desc_ratio = 0.8,
    const float feat_scale_ratio = std::numeric_limits<float>::infinity(),
    const float max_desc_d = std::numeric_limits<float>::infinity()
  ) override
  {
    /*
    const float desc_ratio_2 = desc_ratio*desc_ratio;

    const std::vector<MapLandmark *> & frame_1_3D_pts = frame_1->map_points_;
    const std::vector<MapLandmark *> & frame_2_3D_pts = frame_2->map_points_;
    const size_t frame_1_pts_size = frame_1->getNumberOfFeatures();
    const size_t frame_2_pts_size = frame_2->getNumberOfFeatures();

    features::Regions * const frame_1_regions = frame_1->getRegions();
    features::Regions * const frame_2_regions = frame_2->getRegions();

    bool bCam1_Calibrated = frame_1->getCamCalibrated();
    bool bCam2_Calibrated = frame_2->getCamCalibrated();

    const IntrinsicBase * cam_intrinsic_1 = frame_1->getCameraIntrinsics();
    const IntrinsicBase * cam_intrinsic_2 = frame_2->getCameraIntrinsics();

    // -------------------
    // -- Match unmatched by checking epipolar
    // -------------------
    std::vector<IndexT>putative_matches_1_2_idx(frame_1_pts_size,UndefinedIndexT);

    double startTime = omp_get_wtime();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (IndexT pt_1_i=0; pt_1_i < frame_1_pts_size; ++pt_1_i)
    {
      // If feature is already matched we skip it
      if (frame_1_3D_pts[pt_1_i])
        continue;

      // Get position of point in frame_1
      const Vec2 pt_frame_1 = (bCam1_Calibrated ? frame_1->getFeaturePosition(pt_1_i) : cam_intrinsic_1->remove_disto(frame_1->getFeaturePosition(pt_1_i))).cast<double>();
      // Epipolar line on frame_2
      Vec3 x1_F = F_21 * pt_frame_1.homogeneous() ; // Epipolar line on frame_2

      // Raw descriptors
      void * pt_1_raw_desc;
      void * candidate_pt_desc_raw;

      // Get descriptor from frame_1
      featureExtractor_->getDescriptorRaw(frame_1_regions, pt_1_i, &pt_1_raw_desc);

      //TODO: Get possible candidates through grid
      IndexT best_idx = UndefinedIndexT;
      IndexT second_best_idx = UndefinedIndexT;
      double best_distance_2 = std::numeric_limits<double>::infinity();
      double second_best_distance_2 = std::numeric_limits<double>::infinity();
      double distance_2;

      for (IndexT pt_2_i=0; pt_2_i < frame_2_pts_size; ++pt_2_i)
      {
        // If feature is already matched we skip it
        if (frame_2_3D_pts[pt_2_i])
          continue;

        const Vec2 pt_frame_2 = (bCam2_Calibrated ? frame_2->getFeaturePosition(pt_2_i) : cam_intrinsic_2->remove_disto(frame_2->getFeaturePosition(pt_2_i))).cast<double>();

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_2_regions, pt_2_i, &candidate_pt_desc_raw);

        // Compute distance
        distance_2 = featureExtractor_->SquaredDescriptorDistance(pt_1_raw_desc,candidate_pt_desc_raw);

        // Save if in best two
        if (distance_2 < best_distance_2)
        {
          second_best_distance_2 = best_distance_2;
          second_best_idx = best_idx;
          best_distance_2 = distance_2;
          best_idx = pt_2_i;
        }
        else if (distance_2 < second_best_distance_2)
        {
          second_best_distance_2 = distance_2;
          second_best_idx = pt_2_i;
        }
      }

      // Detect best match
      if (best_idx != UndefinedIndexT && best_distance_2 < max_desc_d_2)
      {
        const Vec2 pt_frame_2 = (bCam2_Calibrated ? frame_2->getFeaturePosition(best_idx) : cam_intrinsic_2->remove_disto(frame_2->getFeaturePosition(best_idx))).cast<double>();

        // Compute distance from epipolar line
        double dF_12 = x1_F.dot(pt_frame_2.homogeneous());
        double d_pt_2_ep_line = (dF_12 * dF_12) /  (x1_F.head<2>().squaredNorm()); // square distance of pt_2 from epipolar line in frame_2

        if (d_pt_2_ep_line > max_epipolar_d2)
          continue;

        if (second_best_idx != UndefinedIndexT  && second_best_distance_2 < max_desc_d_2)
        {
          if ((best_distance_2 / second_best_distance_2) < desc_ratio_2)
          {
            // Best is unique enough
            putative_matches_1_2_idx[pt_1_i] = best_idx;
          }
        }
        else
        {
          // Best is unique
          putative_matches_1_2_idx[pt_1_i] = best_idx;
        }
      }
    }

    std::cout<<"Matching: Epipolar 2D-2D ("<<omp_get_wtime() - startTime<<")\n";

    // -------------------
    // -- Purge matches
    // -------------------
    startTime = omp_get_wtime();

    purgeCandidateMatches(putative_matches_1_2_idx, map_matches_1_2_idx);

    std::cout<<"Matching: Purge ("<<omp_get_wtime() - startTime<<")\n";
*/
  }






};

}
}
