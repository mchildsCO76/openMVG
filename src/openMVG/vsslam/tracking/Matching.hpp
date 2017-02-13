
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once


namespace openMVG  {
namespace VSSLAM  {


  static void purgeCandidateMatches
  (
    std::vector<int> & putative_matches_1_2_idx,
    Hash_Map<size_t,size_t> & matches_1_2_idx
  )
  {
    // Check that two are not matching the same point
    bool bOkMatch = true;
    for (size_t i=0; i<putative_matches_1_2_idx.size(); ++i)
    {
      if (putative_matches_1_2_idx[i] != -1)
      {
        for (size_t j=i+1; j<putative_matches_1_2_idx.size(); ++j)
        {
          // if value is doubled we delete both matches (mathces have to be unique)
          if (putative_matches_1_2_idx[i] == putative_matches_1_2_idx[j])
          {
            bOkMatch = false;
            putative_matches_1_2_idx[j] = -1;
          }
        }
        if (bOkMatch)
        {
          // Match
          matches_1_2_idx[i] = putative_matches_1_2_idx[i];
        }
        else
        {
          // reset duplicate flag
          bOkMatch = true;

        }
      }
    }
  }

  static void purgeCandidateMatches
  (
    const std::vector<MapLandmark *> & frame_1_3D_pts,
    std::vector<int> & putative_matches_1_2_idx,
    Hash_Map<MapLandmark *,size_t> & matches_3D_ptr_cur_idx
  )
  {
    // Check that two are not matching the same point
    bool bOkMatch = true;
    for (size_t i=0; i<putative_matches_1_2_idx.size(); ++i)
    {
      if (putative_matches_1_2_idx[i] != -1)
      {
        for (size_t j=i+1; j<putative_matches_1_2_idx.size(); ++j)
        {
          // if value is doubled we delete both matches (mathces have to be unique)
          if (putative_matches_1_2_idx[i] == putative_matches_1_2_idx[j])
          {
            bOkMatch = false;
            putative_matches_1_2_idx[j] = -1;
          }
        }
        if (bOkMatch)
        {
          // Match
          matches_3D_ptr_cur_idx[frame_1_3D_pts[i]] = putative_matches_1_2_idx[i];
        }
        else
        {
          // reset duplicate flag
          bOkMatch = true;

        }
      }
    }
  }


  // Match two sets of features by:
  //   for each feature in F_1:
  //     compare descriptors only of features in F_2
  //     that are inside a window around position of feature in F_1

  static size_t matching_Window_2D_2D
  (
    Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    Hash_Map<size_t,size_t> & matches_1_2_idx,
    const size_t win_size = 400,
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  )
  {
    features::Regions * frame_1_regions = frame_1->regions_.get();
    const std::vector<Vec2> & frame_1_pts = frame_1->pts_undist_;
    const size_t frame_1_pts_size = frame_1->getNumberOfFeatures();

    features::Regions * frame_2_regions = frame_2->regions_.get();
    const std::vector<Vec2> & frame_2_pts = frame_2->pts_undist_;
    const size_t frame_2_pts_size = frame_2->getNumberOfFeatures();

    // -------------------
    // -- Match reference-current candidates
    // -------------------
    std::vector<int>putative_matches_1_2_idx(frame_1_pts_size,-1);

    double startTime = omp_get_wtime();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t p_i=0; p_i < frame_1_pts_size; ++p_i)
    {
      //TODO: Get possible candidates through grid

      size_t best_idx = std::numeric_limits<size_t>::infinity();
      size_t second_best_idx = std::numeric_limits<size_t>::infinity();
      double best_distance = max_desc_d;//std::numeric_limits<double>::infinity();
      double second_best_distance = max_desc_d;//std::numeric_limits<double>::infinity();

      // Raw descriptors
      void * ref_pt_desc_raw;
      void * candidate_pt_desc_raw;

      // Get descriptor of reference point
      featureExtractor_->getDescriptorRaw(frame_1_regions, p_i, &ref_pt_desc_raw);

      for (size_t c_i=0; c_i < frame_2_pts_size; ++c_i)
      {
        // Check if points are close enough
        if ((frame_2_pts[c_i] - frame_1_pts[p_i]).squaredNorm() > win_size)
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_2_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance
        double distance = featureExtractor_->SquaredDescriptorDistance(ref_pt_desc_raw,candidate_pt_desc_raw);

        // Save if in best two
        if (distance < best_distance)
        {
          second_best_distance = best_distance;
          second_best_idx = best_idx;
          best_distance = distance;
          best_idx = c_i;
        }
        else if (distance < second_best_distance)
        {
          second_best_distance = distance;
          second_best_idx = c_i;
        }
      }

      // Detect best match
      if (best_idx != std::numeric_limits<size_t>::infinity())
      {
        if (second_best_idx != std::numeric_limits<size_t>::infinity())
        {
          if ((best_distance / second_best_distance) < desc_ratio)
          {
            // Best is unique enough
            putative_matches_1_2_idx[p_i] = best_idx;
          }
        }
        else
        {
          // Best is unique
          putative_matches_1_2_idx[p_i] = best_idx;
        }
      }
    }

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Matching: Window 2D-2D ("<<secsElapsed<<")\n";

    // -------------------
    // -- Purge matches
    // -------------------
    startTime = omp_get_wtime();

    purgeCandidateMatches(putative_matches_1_2_idx, matches_1_2_idx);

    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Matching: Purge ("<<secsElapsed<<")\n";

    return matches_1_2_idx.size();
  }

  static size_t matching_AllAll_2D_2D
  (
    Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    Hash_Map<size_t,size_t> & matches_1_2_idx,
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  )
  {

    features::Regions * frame_1_regions = frame_1->regions_.get();
    const size_t frame_1_pts_size = frame_1->getNumberOfFeatures();

    features::Regions * frame_2_regions = frame_2->regions_.get();
    const size_t frame_2_pts_size = frame_2->getNumberOfFeatures();

    // -------------------
    // -- Match reference-current candidates
    // -------------------
    std::vector<int>putative_matches_1_2_idx(frame_1_pts_size,-1);

    double startTime = omp_get_wtime();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t p_i=0; p_i < frame_1_pts_size; ++p_i)
    {
      size_t best_idx = std::numeric_limits<size_t>::infinity();
      size_t second_best_idx = std::numeric_limits<size_t>::infinity();
      double best_distance = max_desc_d;//std::numeric_limits<double>::infinity();
      double second_best_distance = max_desc_d;//std::numeric_limits<double>::infinity();

      // Raw descriptors
      void * ref_pt_desc_raw;
      void * candidate_pt_desc_raw;

      // Get descriptor of reference point
      featureExtractor_->getDescriptorRaw(frame_1_regions, p_i, &ref_pt_desc_raw);

      for (size_t c_i=0; c_i < frame_2_pts_size; ++c_i)
      {
        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_2_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance
        double distance = featureExtractor_->SquaredDescriptorDistance(ref_pt_desc_raw,candidate_pt_desc_raw);

        // Save if in best two
        if (distance < best_distance)
        {
          second_best_distance = best_distance;
          second_best_idx = best_idx;
          best_distance = distance;
          best_idx = c_i;
        }
        else if (distance < second_best_distance)
        {
          second_best_distance = distance;
          second_best_idx = c_i;
        }
      }

      // Detect best match
      if (best_idx != std::numeric_limits<size_t>::infinity())
      {
        if (second_best_idx != std::numeric_limits<size_t>::infinity())
        {
          if ((best_distance / second_best_distance) < desc_ratio)
          {
            // Best is unique enough
            putative_matches_1_2_idx[p_i] = best_idx;
          }
        }
        else
        {
          // Best is unique
          putative_matches_1_2_idx[p_i] = best_idx;
        }
      }
    }

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Matching: AllAll 2D-2D ("<<secsElapsed<<")\n";

    // -------------------
    // -- Purge matches
    // -------------------
    startTime = omp_get_wtime();

    purgeCandidateMatches(putative_matches_1_2_idx, matches_1_2_idx);

    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Matching: Purge ("<<secsElapsed<<")\n";

    return matches_1_2_idx.size();
  }

  // Matching all available 3D points with the 2D features of frame_2
  static size_t matching_AllAll_3D_2D
  (
    Abstract_FeatureExtractor * featureExtractor_,
    const std::vector<MapLandmark *> & frame_1_3D_pts,
    const Frame * frame_2,
    Hash_Map<MapLandmark* ,size_t> &matches_3D_ptr_cur_idx,
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  )
  {
    const size_t frame_1_3D_pts_size = frame_1_3D_pts.size();
    features::Regions * const frame_2_regions = frame_2->regions_.get();
    const size_t frame_2_pts_size = frame_2->getNumberOfFeatures();

    // -------------------
    // -- Match reference-current candidates
    // -------------------
    std::vector<int>matches_1_2_idx(frame_1_3D_pts_size,-1);

    double startTime = omp_get_wtime();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t p_i=0; p_i < frame_1_3D_pts_size; ++p_i)
    {
      MapLandmark * point_3D = frame_1_3D_pts[p_i];
      if (!point_3D)
        continue;

      // Already matched
      if (matches_3D_ptr_cur_idx.find(point_3D) != matches_3D_ptr_cur_idx.end())
        continue;

      size_t best_idx = std::numeric_limits<size_t>::infinity();
      size_t second_best_idx = std::numeric_limits<size_t>::infinity();
      double best_distance = max_desc_d;//std::numeric_limits<double>::infinity();
      double second_best_distance = max_desc_d;//std::numeric_limits<double>::infinity();

      // Raw descriptors
      void * pt_3D_raw_desc = point_3D->bestDesc_;
      void * candidate_pt_desc_raw;

      for (size_t c_i=0; c_i < frame_2_pts_size; ++c_i)
      {
        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_2_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance
        double distance = featureExtractor_->SquaredDescriptorDistance(pt_3D_raw_desc,candidate_pt_desc_raw);

        // Save if in best two
        if (distance < best_distance)
        {
          second_best_distance = best_distance;
          second_best_idx = best_idx;
          best_distance = distance;
          best_idx = c_i;
        }
        else if (distance < second_best_distance)
        {
          second_best_distance = distance;
          second_best_idx = c_i;
        }
      }

      // Detect best match
      if (best_idx != std::numeric_limits<size_t>::infinity())
      {
        if (second_best_idx != std::numeric_limits<size_t>::infinity())
        {
          if ((best_distance / second_best_distance) < desc_ratio)
          {
            // Best is unique enough
            matches_1_2_idx[p_i] = best_idx;
          }
        }
        else
        {
          // Best is unique
          matches_1_2_idx[p_i] = best_idx;
        }
      }
    }

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Matching: AllAll 3D-2D ("<<secsElapsed<<")\n";

    // -------------------
    // -- Purge matches
    // -------------------
    startTime = omp_get_wtime();

    purgeCandidateMatches(frame_1_3D_pts,matches_1_2_idx,matches_3D_ptr_cur_idx);

    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Matching: Purge ("<<secsElapsed<<")\n";

    return matches_3D_ptr_cur_idx.size();

  }

  static size_t matching_Projection_3D_2D
  (
    Abstract_FeatureExtractor * featureExtractor_,
    const std::vector<MapLandmark *> & frame_1_3D_pts,
    const Frame * frame_2,
    Hash_Map<MapLandmark *,size_t> & matches_3D_ptr_cur_idx,
    const size_t win_size = 400, //20*20
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  )
  {
    const size_t frame_1_3D_pts_size = frame_1_3D_pts.size();
    features::Regions * frame_2_regions = frame_2->regions_.get();
    const std::vector<Vec2> & frame_2_pts = frame_2->pts_undist_;

    // -------------------
    // -- Match features of prev_frame(that are already triangulated and in map) with the features of the new frame
    // -------------------
    std::vector<int>matches_1_2_idx(frame_1_3D_pts_size,-1);

    //const Similarity3 & pose_2 = frame_2->pose_;
    const IntrinsicBase * cam_intrinsic_2 = frame_2->getCameraIntrinsics();

    double startTime = omp_get_wtime();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t p_i=0; p_i < frame_1_3D_pts_size; ++p_i)
    {
      MapLandmark * point_3D = frame_1_3D_pts[p_i];

      // Check if its triangulated
      if (!point_3D)
        continue;

      // Already matched
      if (matches_3D_ptr_cur_idx.find(point_3D) != matches_3D_ptr_cur_idx.end())
        continue;

      // Project the point to frame 2
      // TODO: Check how it will work with relative poses
      //Vec3 P_3D_2 = pose_2(point_3D->X_);
      Vec3 P_3D_2;
      getRelativePointPosition(point_3D->X_,point_3D->ref_frame_,P_3D_2,frame_2);
      if (P_3D_2(2) < 0)
        continue;

      // Project the point to image plane of frame 2
      const Vec2 pt_3D_2 = cam_intrinsic_2->cam2ima(cam_intrinsic_2->have_disto()?cam_intrinsic_2->add_disto(P_3D_2.hnormalized()):P_3D_2.hnormalized());

      // Check if projection is actually in the image borders
      if (!frame_2->isPointInFrame(pt_3D_2))
        continue;

      // Raw descriptors
      void * pt_3D_raw_desc = point_3D->bestDesc_;
      void * candidate_pt_desc_raw;

      //TODO: Get possible candidates through grid

      size_t best_idx = std::numeric_limits<size_t>::infinity();
      size_t second_best_idx = std::numeric_limits<size_t>::infinity();
      double best_distance = max_desc_d;//std::numeric_limits<double>::infinity();
      double second_best_distance = max_desc_d;//std::numeric_limits<double>::infinity();

      for (size_t c_i=0; c_i < frame_2_pts.size(); ++c_i)
      {
        // Check if points are close enough
        if ((frame_2_pts[c_i] - pt_3D_2).squaredNorm() > win_size)
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_2_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance
        double distance = featureExtractor_->SquaredDescriptorDistance(pt_3D_raw_desc,candidate_pt_desc_raw);

        // Save if in best two
        if (distance < best_distance)
        {
          second_best_distance = best_distance;
          second_best_idx = best_idx;
          best_distance = distance;
          best_idx = c_i;
        }
        else if (distance < second_best_distance)
        {
          second_best_distance = distance;
          second_best_idx = c_i;
        }
      }

      // Detect best match
      if (best_idx != std::numeric_limits<size_t>::infinity())
      {
        if (second_best_idx != std::numeric_limits<size_t>::infinity())
        {
          if ((best_distance / second_best_distance) < desc_ratio)
          {
            // Best is unique enough
            matches_1_2_idx[p_i] = best_idx;
          }
        }
        else
        {
          // Best is unique
          matches_1_2_idx[p_i] = best_idx;
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

    purgeCandidateMatches(frame_1_3D_pts,matches_1_2_idx,matches_3D_ptr_cur_idx);

    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Matching: Purge ("<<secsElapsed<<")\n";

    return matches_3D_ptr_cur_idx.size();
  }

  static size_t matching_EpipolarLine_2D_2D
  (
    Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    const Mat3 & F_21,
    Hash_Map<size_t,size_t> & matches_1_2_idx,
    const float desc_ratio = 0.8,
    const float max_desc_d = std::numeric_limits<float>::infinity()
  )
  {
    const std::vector<MapLandmark *> & frame_1_3D_pts = frame_1->map_points_;
    const std::vector<MapLandmark *> & frame_2_3D_pts = frame_2->map_points_;
    const size_t frame_1_pts_size = frame_1->getNumberOfFeatures();
    const size_t frame_2_pts_size = frame_2->getNumberOfFeatures();

    features::Regions * frame_1_regions = frame_1->regions_.get();
    features::Regions * frame_2_regions = frame_2->regions_.get();

    bool bCam1_Calibrated = frame_1->getCamCalibrated();
    bool bCam2_Calibrated = frame_2->getCamCalibrated();

    const IntrinsicBase * cam_intrinsic_1 = frame_1->getCameraIntrinsics();
    const IntrinsicBase * cam_intrinsic_2 = frame_2->getCameraIntrinsics();

    // -------------------
    // -- Match unmatched by checking epipolar
    // -------------------

    std::vector<int>putative_matches_1_2_idx(frame_1_pts_size,-1);

    double startTime = omp_get_wtime();

    //#ifdef OPENMVG_USE_OPENMP
    //#pragma omp parallel for schedule(dynamic)
    //#endif
    for (size_t p_i=0; p_i < frame_1_pts_size; ++p_i)
    {
      // If feature is already matched we skip it
      if (frame_1_3D_pts[p_i])
        continue;

      // Get position of point in frame_1
      const Vec2 pt_1 = (bCam1_Calibrated ? frame_1->getFeaturePositionUndistorted(p_i) : cam_intrinsic_1->remove_disto(frame_1->getFeaturePositionUndistorted(p_i))).cast<double>();
      // Epipolar line on frame_2
      Vec3 x1_F = F_21 * pt_1.homogeneous() ; // Epipolar line on frame_2

      // Raw descriptors
      void * pt_1_raw_desc;
      void * candidate_pt_desc_raw;

      // Get descriptor from frame_1
      featureExtractor_->getDescriptorRaw(frame_1_regions, p_i, &pt_1_raw_desc);

      //TODO: Get possible candidates through grid
      size_t best_idx = std::numeric_limits<size_t>::infinity();
      size_t second_best_idx = std::numeric_limits<size_t>::infinity();
      double best_distance = max_desc_d;//std::numeric_limits<double>::infinity();
      double second_best_distance = max_desc_d;//std::numeric_limits<double>::infinity();

      for (size_t c_i=0; c_i < frame_2_pts_size; ++c_i)
      {
        // If feature is already matched we skip it
        if (frame_2_3D_pts[c_i])
          continue;

        const Vec2 pt_2 = (bCam2_Calibrated ? frame_2->getFeaturePositionUndistorted(c_i) : cam_intrinsic_2->remove_disto(frame_2->getFeaturePositionUndistorted(c_i))).cast<double>();

        // Compute distance from epipolar line
        double dF_12 = x1_F.dot(pt_2.homogeneous());
        double d_pt_2_ep_line = (dF_12 * dF_12) /  (x1_F.head<2>().squaredNorm()); // square distance of pt_2 from epipolar line in frame_2

        if (d_pt_2_ep_line > 16)
          continue;

        // Get candidate descriptor
        featureExtractor_->getDescriptorRaw(frame_2_regions, c_i, &candidate_pt_desc_raw);

        // Compute distance
        double distance = featureExtractor_->SquaredDescriptorDistance(pt_1_raw_desc,candidate_pt_desc_raw);

        // Save if in best two
        if (distance < best_distance)
        {
          second_best_distance = best_distance;
          second_best_idx = best_idx;
          best_distance = distance;
          best_idx = c_i;
        }
        else if (distance < second_best_distance)
        {
          second_best_distance = distance;
          second_best_idx = c_i;
        }
      }

      // Detect best match
      if (best_idx != std::numeric_limits<size_t>::infinity())
      {
        if (second_best_idx != std::numeric_limits<size_t>::infinity())
        {
          if ((best_distance / second_best_distance) < desc_ratio)
          {
            // Best is unique enough
            putative_matches_1_2_idx[p_i] = best_idx;
          }
        }
        else
        {
          // Best is unique
          putative_matches_1_2_idx[p_i] = best_idx;
        }
      }
    }

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Matching: Epipolar 2D-2D ("<<secsElapsed<<")\n";

    // -------------------
    // -- Purge matches
    // -------------------
    startTime = omp_get_wtime();

    purgeCandidateMatches(putative_matches_1_2_idx, matches_1_2_idx);

    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Matching: Purge ("<<secsElapsed<<")\n";

    return matches_1_2_idx.size();

  }



}
}
