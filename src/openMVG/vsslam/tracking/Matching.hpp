
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once


namespace openMVG  {
namespace VSSLAM  {

  // Match two sets of features by:
  //   for each feature in F_1:
  //     compare descriptors only of features in F_2
  //     that are inside a window around position of feature in F_1

  static size_t matchingByWindow_2D_2D
  (
    Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    Hash_Map<size_t,size_t> & match_1_2_idx,
    const size_t win_size = 50,
    const float desc_ratio = 0.8
  )
  {
    features::Regions * frame_1_regions = frame_1->regions_.get();
    const std::vector<Vec2> & frame_1_pts = frame_1->pts_undist_;

    features::Regions * frame_2_regions = frame_2->regions_.get();
    const std::vector<Vec2> & frame_2_pts = frame_2->pts_undist_;

    // -------------------
    // -- Match reference-current candidates
    // -------------------
    std::cout<<"Matching by window search \n";

    std::vector<int>matches_1_2_idx(frame_1_pts.size(),-1);

    double startTime = omp_get_wtime();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t p_i=0; p_i < frame_1_pts.size(); ++p_i)
    {
      //TODO: Get possible candidates through grid

      size_t best_idx = std::numeric_limits<size_t>::infinity();
      size_t second_best_idx = std::numeric_limits<size_t>::infinity();
      double best_distance = 30;//std::numeric_limits<double>::infinity();
      double second_best_distance = 30;//std::numeric_limits<double>::infinity();

      // Raw descriptors
      void * ref_pt_desc_raw;
      void * candidate_pt_desc_raw;

      // Get descriptor of reference point
      featureExtractor_->getDescriptorRaw(frame_1_regions, p_i, &ref_pt_desc_raw);

      for (size_t c_i=0; c_i < frame_2_pts.size(); ++c_i)
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
    std::cout<<"Matching by window search:"<<secsElapsed<<"\n";


    // -------------------
    // -- Prune matches and insert
    // -------------------
    std::cout<<"Purging candidates and copy\n";

    startTime = omp_get_wtime();

    // Check that two are not matching the same point
    bool bOkMatch = true;
    for (size_t i=0; i<matches_1_2_idx.size(); ++i)
    {
      if (matches_1_2_idx[i] != -1)
      {
        for (size_t j=i+1; j<matches_1_2_idx.size(); ++j)
        {
          // if value is doubled we delete both matches (mathces have to be unique)
          if (matches_1_2_idx[i] == matches_1_2_idx[j])
          {
            bOkMatch = false;
            matches_1_2_idx[j] = -1;
          }
        }
        if (bOkMatch)
        {
          // Match
          match_1_2_idx[i] = matches_1_2_idx[i];
        }
        else
        {
          // reset duplicate flag
          bOkMatch = true;

        }
      }
    }

    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Purging time:"<<secsElapsed<<"\n";

    return match_1_2_idx.size();
  }


  static size_t matching_AllAll_2D_2D
  (
    Abstract_FeatureExtractor * featureExtractor_,
    const Frame * frame_1,
    const Frame * frame_2,
    Hash_Map<size_t,size_t> & match_1_2_idx,
    const float desc_ratio = 0.8
  )
  {

    features::Regions * frame_1_regions = frame_1->regions_.get();
    const std::vector<Vec2> & frame_1_pts = frame_1->pts_undist_;
    const size_t frame_1_pts_size = frame_1_pts.size();

    features::Regions * frame_2_regions = frame_2->regions_.get();
    const std::vector<Vec2> & frame_2_pts = frame_2->pts_undist_;

    // -------------------
    // -- Match reference-current candidates
    // -------------------
    std::cout<<"Matching AllAll_2D_2D \n";

    std::vector<int>matches_1_2_idx(frame_1_pts_size,-1);

    double startTime = omp_get_wtime();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t p_i=0; p_i < frame_1_pts_size; ++p_i)
    {
      size_t best_idx = std::numeric_limits<size_t>::infinity();
      size_t second_best_idx = std::numeric_limits<size_t>::infinity();
      double best_distance = 30;//std::numeric_limits<double>::infinity();
      double second_best_distance = 30;//std::numeric_limits<double>::infinity();

      // Raw descriptors
      void * ref_pt_desc_raw;
      void * candidate_pt_desc_raw;

      // Get descriptor of reference point
      featureExtractor_->getDescriptorRaw(frame_1_regions, p_i, &ref_pt_desc_raw);

      for (size_t c_i=0; c_i < frame_2_pts.size(); ++c_i)
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
    std::cout<<"Matching AllAll_2D_2D:"<<secsElapsed<<"\n";


    // -------------------
    // -- Prune matches and insert
    // -------------------
    std::cout<<"Purging candidates and copy\n";

    startTime = omp_get_wtime();

    // Check that two are not matching the same point
    bool bOkMatch = true;
    for (size_t i=0; i<matches_1_2_idx.size(); ++i)
    {
      if (matches_1_2_idx[i] != -1)
      {
        for (size_t j=i+1; j<matches_1_2_idx.size(); ++j)
        {
          // if value is doubled we delete both matches (mathces have to be unique)
          if (matches_1_2_idx[i] == matches_1_2_idx[j])
          {
            bOkMatch = false;
            matches_1_2_idx[j] = -1;
          }
        }
        if (bOkMatch)
        {
          // Match
          match_1_2_idx[i] = matches_1_2_idx[i];
        }
        else
        {
          // reset duplicate flag
          bOkMatch = true;

        }
      }
    }

    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Purging time:"<<secsElapsed<<"\n";

    return match_1_2_idx.size();
  }




  static size_t matchingAllAll_3D_2D
  (
    Abstract_FeatureExtractor * featureExtractor_,
    const std::vector<MapLandmark *> & frame_1_3D_pts,
    const Frame * frame_2,
    Hash_Map<MapLandmark* ,size_t> &matches_3D_ptr_cur_idx,
    const float desc_ratio = 0.8
  )
  {
    const size_t frame_1_3D_pts_size = frame_1_3D_pts.size();
    features::Regions * frame_2_regions = frame_2->regions_.get();
    const std::vector<Vec2> & frame_2_pts = frame_2->pts_undist_;

    // -------------------
    // -- Match reference-current candidates
    // -------------------
    std::cout<<"Matching All-All search \n";

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
      double best_distance = 30;//std::numeric_limits<double>::infinity();
      double second_best_distance = 30;//std::numeric_limits<double>::infinity();

      // Raw descriptors
      void * pt_3D_raw_desc = point_3D->bestDesc_;
      void * candidate_pt_desc_raw;

      for (size_t c_i=0; c_i < frame_2_pts.size(); ++c_i)
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
    std::cout<<"Matching All-All search:"<<secsElapsed<<"\n";


    // -------------------
    // -- Prune matches and insert
    // -------------------
    std::cout<<"Purging candidates and copy\n";

    startTime = omp_get_wtime();

    // Check that two are not matching the same point
    bool bOkMatch = true;
    for (size_t i=0; i<matches_1_2_idx.size(); ++i)
    {
      if (matches_1_2_idx[i] != -1)
      {
        for (size_t j=i+1; j<matches_1_2_idx.size(); ++j)
        {
          // if value is doubled we delete both matches (mathces have to be unique)
          if (matches_1_2_idx[i] == matches_1_2_idx[j])
          {
            bOkMatch = false;
            matches_1_2_idx[j] = -1;
          }
        }
        if (bOkMatch)
        {
          // Match
          matches_3D_ptr_cur_idx[frame_1_3D_pts[i]] = matches_1_2_idx[i];
        }
        else
        {
          // reset duplicate flag
          bOkMatch = true;

        }
      }
    }

    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Purging time:"<<secsElapsed<<"\n";

    return matches_3D_ptr_cur_idx.size();

  }

  static size_t matchingByProjection
  (
    Abstract_FeatureExtractor * featureExtractor_,
    const std::vector<MapLandmark *> & frame_1_3D_pts,
    const Frame * frame_2,
    Hash_Map<MapLandmark *,size_t> & matches_3D_ptr_cur_idx,
    const size_t win_size = 50,
    const float desc_ratio = 0.8
  )
  {
    features::Regions * frame_2_regions = frame_2->regions_.get();
    const std::vector<Vec2> & frame_2_pts = frame_2->pts_undist_;

    // -------------------
    // -- Match features of prev_frame(that are already triangulated and in map) with the features of the new frame
    // -------------------
    std::cout<<"Matching by projection search \n";

    std::vector<int>matches_1_2_idx(frame_1_3D_pts.size(),-1);


    const Similarity3 & pose_2 = frame_2->pose_;
    const Camera * cam_1 = frame_2->cam_;
    const IntrinsicBase * cam_intrinsic_1 = cam_1->cam_intrinsic_ptr;
    const Camera * cam_2 = frame_2->cam_;
    const IntrinsicBase * cam_intrinsic_2 = cam_2->cam_intrinsic_ptr;

    std::cout<<"AA\n";
    double startTime = omp_get_wtime();

    #ifdef OPENMVG_USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (size_t p_i=0; p_i < frame_1_3D_pts.size(); ++p_i)
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
      Vec3 P_3D_2 = pose_2(point_3D->X_);

      if (P_3D_2[2] < 0)
        continue;

      // Project the point to image plane of frame 2
      const Vec2 pt_3D_2 = cam_intrinsic_2->cam2ima(cam_intrinsic_2->have_disto()?cam_intrinsic_2->add_disto(P_3D_2.hnormalized()):P_3D_2.hnormalized());

      // Raw descriptors
      void * pt_3D_raw_desc = point_3D->bestDesc_;
      void * candidate_pt_desc_raw;

      //TODO: Get possible candidates through grid

      size_t best_idx = std::numeric_limits<size_t>::infinity();
      size_t second_best_idx = std::numeric_limits<size_t>::infinity();
      double best_distance = 30;//std::numeric_limits<double>::infinity();
      double second_best_distance = 30;//std::numeric_limits<double>::infinity();

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
    std::cout<<"Matching by projection search:"<<secsElapsed<<"\n";


    // -------------------
    // -- Prune matches and insert
    // -------------------
    std::cout<<"Purging candidates and copy\n";

    startTime = omp_get_wtime();

    // Check that two are not matching the same point
    bool bOkMatch = true;
    for (size_t i=0; i<matches_1_2_idx.size(); ++i)
    {
      if (matches_1_2_idx[i] != -1)
      {
        for (size_t j=i+1; j<matches_1_2_idx.size(); ++j)
        {
          // if value is doubled we delete both matches (mathces have to be unique)
          if (matches_1_2_idx[i] == matches_1_2_idx[j])
          {
            bOkMatch = false;
            matches_1_2_idx[j] = -1;
          }
        }
        if (bOkMatch)
        {
          // Match
          matches_3D_ptr_cur_idx[frame_1_3D_pts[i]] = matches_1_2_idx[i];
        }
        else
        {
          // reset duplicate flag
          bOkMatch = true;

        }
      }
    }

    stopTime = omp_get_wtime();
    secsElapsed = stopTime - startTime; // that's all !
    std::cout<<"Purging time:"<<secsElapsed<<"\n";

    return matches_3D_ptr_cur_idx.size();
  }

}
}
