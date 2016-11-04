
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/sfm/pipelines/sfm_engine.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/tracks/tracks.hpp"

#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/histogram/histogram.hpp"

namespace openMVG {
namespace sfm {

struct SlamPP_Data
{
  IndexT next_free_id_slamPP = 0;

  /// SlamPP index mapping
  Hash_Map<IndexT,IndexT> camera_ids_omvg_slamPP;
  Hash_Map<IndexT,IndexT> track_ids_omvg_slamPP;
  Hash_Map<IndexT,IndexT> camera_ids_slamPP_omvg;
  Hash_Map<IndexT,IndexT> track_ids_slamPP_omvg;

  std::string slamPP_dataset_filename;
  std::ofstream slamPP_DatasetFile;
  
  void createLogFile()
  {
    slamPP_DatasetFile.open( slamPP_dataset_filename.c_str(), std::ios::out );
  }
  void closeLogFile()
  {
    slamPP_DatasetFile.flush();
    slamPP_DatasetFile.close();
  }

  // Get Index Of Camera
  bool getCamId_SlamPP(IndexT &camId_omvg, IndexT &camId_slamPP)
  {
    // Check if slamPP camera id exists
    if (camera_ids_omvg_slamPP.count(camId_omvg) == 0)
      return false;
    
    camId_slamPP = camera_ids_omvg_slamPP[camId_omvg];
    return true;
  }
  bool getCamId_OpenMVG(IndexT &camId_omvg, IndexT &camId_slamPP)
  {
    // Check if omvg camera id exists
    if (camera_ids_slamPP_omvg.count(camId_slamPP) == 0)
      return false;
    
    camId_omvg = camera_ids_slamPP_omvg[camId_slamPP];
    return true;
  }
  // Set Index of Camera
  void setCamId_SlamPP(const IndexT &camId_omvg, const IndexT &camId_slamPP)
  {    
    camera_ids_omvg_slamPP[camId_omvg] = camId_slamPP;
    camera_ids_slamPP_omvg[camId_slamPP] = camId_omvg;
  }

  
  bool getTrackId_SlamPP(IndexT &trackId_omvg, IndexT &trackId_slamPP)
  {
    // Check if slampp track id exists
    if (track_ids_omvg_slamPP.count(trackId_omvg) == 0)
      return false;
    
    trackId_slamPP = track_ids_omvg_slamPP[trackId_omvg];
    return true;
  }
  bool getTrackId_OpenMVG(IndexT &trackId_omvg, IndexT &trackId_slamPP)
  {
    // Check if omvg track id exists
    if (track_ids_slamPP_omvg.count(trackId_slamPP) == 0)
      return false;
    
    trackId_omvg = track_ids_slamPP_omvg[trackId_slamPP];
    return true;
  }
  void setTrackId_SlamPP(const IndexT &trackId_omvg, const IndexT &trackId_slamPP)
  {    
    track_ids_omvg_slamPP[trackId_omvg] = trackId_slamPP;
    track_ids_slamPP_omvg[trackId_slamPP] = trackId_omvg;
  }
  
  IndexT getNextFreeSlamPPId()
  {
    IndexT freeId = next_free_id_slamPP;
    next_free_id_slamPP++;
    return freeId;
  }
  
};

/// Incremental SfM Pipeline Reconstruction Engine.
class IncrementalSfMReconstructionEngine : public ReconstructionEngine
{
public:

  IncrementalSfMReconstructionEngine(
    const SfM_Data & sfm_data,
    const std::string & soutDirectory,
    const std::string & loggingFile = "");

  ~IncrementalSfMReconstructionEngine();

  void SetFeaturesProvider(Features_Provider * provider);
  void SetMatchesProvider(Matches_Provider * provider);

  virtual bool Process();

  void setInitialPair(const Pair & initialPair)
  {
    initial_pair_ = initialPair;
  }

  /// Initialize tracks
  bool InitLandmarkTracks();

  /// Select a candidate initial pair
  bool ChooseInitialPair(Pair & initialPairIndex) const;

  /// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
  bool MakeInitialPair3D(const Pair & initialPair);

  /// Automatic initial pair selection (based on a 'baseline' computation score)
  bool AutomaticInitialPairChoice(Pair & initialPair) const;

  /**
   * Set the default lens distortion type to use if it is declared unknown
   * in the intrinsics camera parameters by the previous steps.
   *
   * It can be declared unknown if the type cannot be deduced from the metadata.
   */
  void SetUnknownCameraType(const cameras::EINTRINSIC camType)
  {
    cam_type_ = camType;
  }

  /// Set which BA are performed
  void setBAOptimizations
  (
    bool performGlobalBA = false,
    bool performInitialTwoViewBA = true,
    bool performLocalPoseBA = true
  )
  {
    performGlobalBA_ = performGlobalBA;
    performInitialTwoViewBA_ = performInitialTwoViewBA;
    performLocalPoseBA_ = performLocalPoseBA;
  }
  /// Set if global outlier removal procedure is done
  void setGlobalOutlierRemoval(bool &performGOR)
  {
    performGlobalOutlierRemoval_ = performGOR;
  }
  /// Set if consistency check of incremental log is done
  void setConsistencyCheck(bool &performCC)
  {
    performConsistencyCheck_ = performCC;
  }

  void setSlamPPOutputFile(std::string &filename)
  {
    slam_pp_data.slamPP_dataset_filename = filename;
  }

protected:


private:

  /// Return MSE (Mean Square Error) and a histogram of residual values.
  double ComputeResidualsHistogram(Histogram<double> * histo);

  /// List the images that the greatest number of matches to the current 3D reconstruction.
  bool FindImagesWithPossibleResection(std::vector<size_t> & vec_possible_indexes);

  /// Add a single Image to the scene and triangulate new possible tracks.
  bool Resection(const size_t imageIndex);

  /// Bundle adjustment to refine Structure; Motion and Intrinsics
  bool BundleAdjustment();

  /// Discard track with too large residual error
  bool badTrackRejector(double dPrecision, size_t count = 0);

  void resetIncrementStep();

  bool checkIncrementalConsistency();

  void ExportIncrementToTextFile_SlamPP();

  //----
  //-- Data
  //----

  // HTML logger
  std::shared_ptr<htmlDocument::htmlDocumentStream> html_doc_stream_;
  std::string sLogging_file_;

  // Parameter
  Pair initial_pair_;
  cameras::EINTRINSIC cam_type_; // The camera type for the unknown cameras

  //-- Data provider
  Features_Provider  * features_provider_;
  Matches_Provider  * matches_provider_;

  // Incremental change
  std::pair<std::set<IndexT>, std::set<IndexT> > increment_camera_;
  std::pair<std::set<IndexT>, std::set<IndexT> > increment_structure_;
  std::pair< Hash_Map<IndexT, std::set<IndexT> >, Hash_Map<IndexT, std::set<IndexT> > > increment_observation_;
  // History
  std::vector<std::pair<std::set<IndexT>, std::set<IndexT> > > history_i_camera_;
  std::vector<std::pair<std::set<IndexT>, std::set<IndexT> > > history_i_structure_;
  std::vector<std::pair< Hash_Map<IndexT, std::set<IndexT> >, Hash_Map<IndexT, std::set<IndexT> > > > history_i_observations_;

  // Owner of tracks
  Hash_Map<IndexT,IndexT> owner_track_cam_id;

  // Temporary data
  openMVG::tracks::STLMAPTracks map_tracks_; // putative landmark tracks (visibility per 3D point)
  Hash_Map<IndexT, double> map_ACThreshold_; // Per camera confidence (A contrario estimated threshold error)

  std::set<size_t> set_remaining_view_id_;     // Remaining camera index that can be used for resection

  // SlamPP test data properties
  float dBAThresholdGroup_;                   // Insert group of images to BA together (if more than dBAThresholdGroup% of points that the best candidate has)

  bool performGlobalBA_;
  bool performInitialTwoViewBA_;
  bool performLocalPoseBA_;
  bool performGlobalOutlierRemoval_;

  bool performConsistencyCheck_;
  
  // SLAM++
  SlamPP_Data slam_pp_data;
};

} // namespace sfm
} // namespace openMVG

