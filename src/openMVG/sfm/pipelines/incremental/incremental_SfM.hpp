
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
#include "openMVG/sfm/pipelines/incremental/SlamPP_GraphFileExport.hpp"

namespace openMVG {
namespace sfm {


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
    const bool performGlobalBA = false,
    const bool performInitialTwoViewBA = true,
    const bool performLocalPoseBA = true
  )
  {
    performGlobalBA_ = performGlobalBA;
    performInitialTwoViewBA_ = performInitialTwoViewBA;
    performLocalPoseBA_ = performLocalPoseBA;
  }
  /// Set if global outlier removal procedure is done
  void setOutlierRemoval(const bool &performGOR, const bool &performLOR)
  {
    performGlobalOutlierRemoval_ = performGOR;
    performLocalOutlierRemoval_ = performLOR;
    
  std::cout<<"OUTLIUE C: "<<performGlobalOutlierRemoval_<<"\n";
  }
  /// Set if consistency check of incremental log is done
  void setConsistencyCheck(const bool &performCC)
  {
    performConsistencyCheck_ = performCC;
  }

  void setGraphFileOutputFile(const std::string &filename)
  {
    slam_pp_data.slamPP_dataset_filename = filename;
  }
  void setGraphVertexOutputTypes(const int camVertexOutputType, const int landmarkOutputType)
  {
    slam_pp_data.iOutputVertexType = camVertexOutputType; // 0 - SE3; 1 - Sim3
    slam_pp_data.iOutputLandmarkType = landmarkOutputType; // 0 - Eucliean (world); 1 - inverse depth (reference)
  }

  void setExportTwoFoldGraphFile(const bool bTF_Graph)
  {
    slam_pp_data.bTwoFoldGraphFile = bTF_Graph;
  }
  
  void ExportIncrementToGraphFile_SlamPP();
  void ExportTwoFoldIncrementToGraphFile_SlamPP();
  void ExportTwoFoldTotalProcessByIncrementToGraphFile_SlamPP();
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

  // Data structure for removed tracks and measurements
  Hash_Map<IndexT, IndexT> structure_last_removed;
  Hash_Map<IndexT, Hash_Map<IndexT, IndexT> > observation_last_removed;


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
  bool performLocalOutlierRemoval_;

  bool performConsistencyCheck_;
  
  // SLAM++
  SlamPP_Data slam_pp_data;
  IndexT current_recon_iteration;
};

} // namespace sfm
} // namespace openMVG

