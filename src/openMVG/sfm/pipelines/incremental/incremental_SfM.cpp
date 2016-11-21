
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"
#include "openMVG/sfm/pipelines/incremental/incremental_SfM.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/cameras/cameras.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"
#include "openMVG/sfm/pipelines/localization/SfM_Localizer.hpp"

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"
#include "openMVG/graph/connectedComponent.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/htmlDoc/htmlDoc.hpp"
#include "third_party/progress/progress.hpp"
#include "ceres/rotation.h"
#include <lemon/bfs.h>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

namespace openMVG {
namespace sfm {

using namespace openMVG::geometry;
using namespace openMVG::cameras;

IncrementalSfMReconstructionEngine::IncrementalSfMReconstructionEngine(
  const SfM_Data & sfm_data,
  const std::string & soutDirectory,
  const std::string & sloggingFile)
  : ReconstructionEngine(sfm_data, soutDirectory),
    sLogging_file_(sloggingFile),
    initial_pair_(Pair(0,0)),
    cam_type_(EINTRINSIC(PINHOLE_CAMERA_RADIAL3))
{
  if (!sLogging_file_.empty())
  {
    // setup HTML logger
    html_doc_stream_ = std::make_shared<htmlDocument::htmlDocumentStream>("SequentialReconstructionEngine SFM report.");
    html_doc_stream_->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("SequentialSfMReconstructionEngine")));
    html_doc_stream_->pushInfo("<hr>");

    html_doc_stream_->pushInfo( "Dataset info:");
    html_doc_stream_->pushInfo( "Views count: " +
      htmlDocument::toString( sfm_data.GetViews().size()) + "<br>");
  }
  // Init remaining image list
  for (Views::const_iterator itV = sfm_data.GetViews().begin();
    itV != sfm_data.GetViews().end(); ++itV)
  {
    set_remaining_view_id_.insert(itV->second.get()->id_view);
  }

  dBAThresholdGroup_ = 0.5f;

  performGlobalBA_ = false;
  performInitialTwoViewBA_ = true;
  performLocalPoseBA_ = true;
  performGlobalOutlierRemoval_ = false;
  performConsistencyCheck_ = true;

  performLocalOutlierRemoval_ = false;


  bOrderedProcessing_ = true;
  orderWindowSize_ = 1;
  bTryAllViews_ = true;

}

IncrementalSfMReconstructionEngine::~IncrementalSfMReconstructionEngine()
{
  if (!sLogging_file_.empty())
  {
    // Save the reconstruction Log
    std::ofstream htmlFileStream(sLogging_file_.c_str());
    htmlFileStream << html_doc_stream_->getDoc();
  }
}

void IncrementalSfMReconstructionEngine::SetFeaturesProvider(Features_Provider * provider)
{
  features_provider_ = provider;
}

void IncrementalSfMReconstructionEngine::SetMatchesProvider(Matches_Provider * provider)
{
  matches_provider_ = provider;
}

bool IncrementalSfMReconstructionEngine::Process() {

  //-------------------
  //-- Incremental reconstruction
  //-------------------
  current_recon_iteration = 0;

  if (!InitLandmarkTracks())
    return false;

  // Initial pair choice
  if (initial_pair_ == Pair(0,0))
  {
    if (!AutomaticInitialPairChoice(initial_pair_))
    {
      // Cannot find a valid initial pair, try to set it by hand?
      if (!ChooseInitialPair(initial_pair_))
      {
        return false;
      }
    }
  }
  // Else a starting pair was already initialized before

  // Initial pair Essential Matrix and [R|t] estimation.
  if (!MakeInitialPair3D(initial_pair_))
    return false;

  // Initial Incremental step
  current_recon_iteration++;
  resetIncrementStep();

  if (performConsistencyCheck_)
    checkIncrementalConsistency();

  // Compute robust Resection of remaining images
  // - group of images will be selected and resection + scene completion will be tried
  size_t resectionGroupIndex = 0;
  std::vector<size_t> vec_possible_resection_indexes;

// bOrderedProcessing_
  while (FindImagesWithPossibleResection(vec_possible_resection_indexes))
  {
    std::cout<<"Found views: "<<vec_possible_resection_indexes.size()<<"\n";
    bool bImageAdded = false;
    // Add images to the 3D reconstruction
    for (std::vector<size_t>::const_iterator iter = vec_possible_resection_indexes.begin();
      iter != vec_possible_resection_indexes.end(); ++iter)
    {
      std::cout<<"Try resection: "<<*iter<<"\n";
      bool bCurrentImageAdded = Resection(*iter);
      bImageAdded |= bCurrentImageAdded;
      if (bCurrentImageAdded)
      {
        set_remaining_view_id_.erase(*iter);
        // After every image its a new iteration
        current_recon_iteration++;
        resetIncrementStep();
        std::cout<<"remove view: "<<*iter<<"\n";
      }
      if (*iter > lastUsedViewId_)
        lastUsedViewId_ = *iter;
    }

    if (bImageAdded)
    {
      std::cout<<"Scene added \n";
      // Scene logging as ply for visual debug
      std::ostringstream os;
      os << std::setw(8) << std::setfill('0') << resectionGroupIndex << "_Resection";
      Save(sfm_data_, stlplus::create_filespec(sOut_directory_, os.str(), ".ply"), ESfM_Data(ALL));

      // Perform BA until all point are under the given precision
      if (performGlobalBA_)
      {
        do
        {
          BundleAdjustment();
        }
        while (performGlobalOutlierRemoval_ && badTrackRejector(16.0, 50));
        //if (performGlobalOutlierRemoval_)
          //eraseUnstablePosesAndObservations(sfm_data_);
      }
      else if(performGlobalOutlierRemoval_)
      {
        badTrackRejector(16.0, 50);
      }
      else
      {
        // Perform just most wide outlier removal to remove far far away points
        badTrackRejector(64.0, 0);
      }

      if (performConsistencyCheck_)
        checkIncrementalConsistency();

    }
    ++resectionGroupIndex;
  }

/*
  // Ensure there is no remaining outliers
  if (performGlobalOutlierRemoval_)
  {
    if (badTrackRejector(4.0, 0))
    {
      //eraseUnstablePosesAndObservations(sfm_data_);
    }
  }
*/

  // Go through incremental logs and put everything to file
  ExportTwoFoldTotalProcessByIncrementToGraphFile_SlamPP();


  //-- Reconstruction done.

  //-- Display some statistics
  std::cout << "\n\n-------------------------------" << "\n"
    << "-- Structure from Motion (statistics):\n"
    << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
    << " from " << sfm_data_.GetViews().size() << " input images.\n"
    << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "\n"
    << "-------------------------------" << "\n";

  Histogram<double> h;
  ComputeResidualsHistogram(&h);
  std::cout << "\nHistogram of residuals:" << h.ToString() << std::endl;

  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Structure from Motion process finished.";
    html_doc_stream_->pushInfo("<hr>");
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- Structure from Motion (statistics):<br>"
      << "-- #Camera calibrated: " << sfm_data_.GetPoses().size()
      << " from " <<sfm_data_.GetViews().size() << " input images.<br>"
      << "-- #Tracks, #3D points: " << sfm_data_.GetLandmarks().size() << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());

    html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of reprojection-residuals"));

    const std::vector<double> xBin = h.GetXbinsValue();
    std::pair< std::pair<double,double>, std::pair<double,double> > range =
      autoJSXGraphViewport<double>(xBin, h.GetHist());

    htmlDocument::JSXGraphWrapper jsxGraph;
    jsxGraph.init("3DtoImageResiduals",600,300);
    jsxGraph.addXYChart(xBin, h.GetHist(), "line,point");
    jsxGraph.UnsuspendUpdate();
    jsxGraph.setViewport(range);
    jsxGraph.close();
    html_doc_stream_->pushInfo(jsxGraph.toStr());
  }
  return true;
}

/// Select a candidate initial pair
bool IncrementalSfMReconstructionEngine::ChooseInitialPair(Pair & initialPairIndex) const
{
  if (initial_pair_ != Pair(0,0))
  {
    // Internal initial pair is already initialized (so return it)
    initialPairIndex = initial_pair_;
  }
  else
  {
    // List Views that supports valid intrinsic
    std::set<IndexT> valid_views;
    for (Views::const_iterator it = sfm_data_.GetViews().begin();
      it != sfm_data_.GetViews().end(); ++it)
    {
      const View * v = it->second.get();
      if( sfm_data_.GetIntrinsics().find(v->id_intrinsic) != sfm_data_.GetIntrinsics().end())
        valid_views.insert(v->id_view);
    }

    if (sfm_data_.GetIntrinsics().empty() || valid_views.empty())
    {
      std::cerr
        << "There is no defined intrinsic data in order to compute an essential matrix for the initial pair."
        << std::endl;
      return false;
    }

    std::cout << std::endl
      << "----------------------------------------------------\n"
      << "SequentialSfMReconstructionEngine::ChooseInitialPair\n"
      << "----------------------------------------------------\n"
      << " Pairs that have valid intrinsic and high support of points are displayed:\n"
      << " Choose one pair manually by typing the two integer indexes\n"
      << "----------------------------------------------------\n"
      << std::endl;

    // Try to list the 10 top pairs that have:
    //  - valid intrinsics,
    //  - valid estimated Fundamental matrix.
    std::vector< size_t > vec_NbMatchesPerPair;
    std::vector<openMVG::matching::PairWiseMatches::const_iterator> vec_MatchesIterator;
    const openMVG::matching::PairWiseMatches & map_Matches = matches_provider_->pairWise_matches_;
    for (openMVG::matching::PairWiseMatches::const_iterator
      iter = map_Matches.begin();
      iter != map_Matches.end(); ++iter)
    {
      const Pair current_pair = iter->first;
      if (valid_views.count(current_pair.first) &&
        valid_views.count(current_pair.second) )
      {
        vec_NbMatchesPerPair.push_back(iter->second.size());
        vec_MatchesIterator.push_back(iter);
      }
    }
    // sort the Pairs in descending order according their correspondences count
    using namespace stl::indexed_sort;
    std::vector< sort_index_packet_descend< size_t, size_t> > packet_vec(vec_NbMatchesPerPair.size());
    sort_index_helper(packet_vec, &vec_NbMatchesPerPair[0], std::min((size_t)10, vec_NbMatchesPerPair.size()));

    for (size_t i = 0; i < std::min((size_t)10, vec_NbMatchesPerPair.size()); ++i) {
      const size_t index = packet_vec[i].index;
      openMVG::matching::PairWiseMatches::const_iterator iter = vec_MatchesIterator[index];
      std::cout << "(" << iter->first.first << "," << iter->first.second <<")\t\t"
        << iter->second.size() << " matches" << std::endl;
    }

    // Ask the user to choose an initial pair (by set some view ids)
    std::cout << std::endl << " type INITIAL pair ids: X enter Y enter\n";
    int val, val2;
    if ( std::cin >> val && std::cin >> val2) {
      initialPairIndex.first = val;
      initialPairIndex.second = val2;
    }
  }

  std::cout << "\nPutative starting pair is: (" << initialPairIndex.first
      << "," << initialPairIndex.second << ")" << std::endl;

  // Check validity of the initial pair indices:
  if (features_provider_->feats_per_view.find(initialPairIndex.first) == features_provider_->feats_per_view.end() ||
      features_provider_->feats_per_view.find(initialPairIndex.second) == features_provider_->feats_per_view.end())
  {
    std::cerr << "At least one of the initial pair indices is invalid."
      << std::endl;
    return false;
  }
  return true;
}

bool IncrementalSfMReconstructionEngine::InitLandmarkTracks()
{
  // Compute tracks from matches
  tracks::TracksBuilder tracksBuilder;

  {
    // List of features matches for each couple of images
    const openMVG::matching::PairWiseMatches & map_Matches = matches_provider_->pairWise_matches_;
    std::cout << "\n" << "Track building" << std::endl;

    tracksBuilder.Build(map_Matches);
    std::cout << "\n" << "Track filtering" << std::endl;
    tracksBuilder.Filter();
    std::cout << "\n" << "Track export to internal struct" << std::endl;
    //-- Build tracks with STL compliant type :
    tracksBuilder.ExportToSTL(map_tracks_);

    std::cout << "\n" << "Track stats" << std::endl;
    {
      std::ostringstream osTrack;
      //-- Display stats :
      //    - number of images
      //    - number of tracks
      std::set<size_t> set_imagesId;
      tracks::TracksUtilsMap::ImageIdInTracks(map_tracks_, set_imagesId);
      osTrack << "------------------" << "\n"
        << "-- Tracks Stats --" << "\n"
        << " Tracks number: " << tracksBuilder.NbTracks() << "\n"
        << " Images Id: " << "\n";
      std::copy(set_imagesId.begin(),
        set_imagesId.end(),
        std::ostream_iterator<size_t>(osTrack, ", "));
      osTrack << "\n------------------" << "\n";

      std::map<size_t, size_t> map_Occurence_TrackLength;
      tracks::TracksUtilsMap::TracksLength(map_tracks_, map_Occurence_TrackLength);
      osTrack << "TrackLength, Occurrence" << "\n";
      for (std::map<size_t, size_t>::const_iterator iter = map_Occurence_TrackLength.begin();
        iter != map_Occurence_TrackLength.end(); ++iter)  {
        osTrack << "\t" << iter->first << "\t" << iter->second << "\n";
      }
      osTrack << "\n";
      std::cout << osTrack.str();
    }
  }
  return map_tracks_.size() > 0;
}

bool IncrementalSfMReconstructionEngine::AutomaticInitialPairChoice(Pair & initial_pair) const
{
  // select a pair that have the largest baseline (mean angle between it's bearing vectors).

  const unsigned iMin_inliers_count = 100;
  const float fRequired_min_angle = 3.0f;
  const float fLimit_max_angle = 60.0f; // More than 60 degree, we cannot rely on matches for initial pair seeding

  // List Views that support valid intrinsic (view that could be used for Essential matrix computation)
  std::set<IndexT> valid_views;
  for (Views::const_iterator it = sfm_data_.GetViews().begin();
    it != sfm_data_.GetViews().end(); ++it)
  {
    const View * v = it->second.get();
    if (sfm_data_.GetIntrinsics().count(v->id_intrinsic))
      valid_views.insert(v->id_view);
  }

  if (valid_views.size() < 2)
  {
    return false; // There is not view that support valid intrinsic data
  }

  std::vector<std::pair<double, Pair> > scoring_per_pair;

  // Compute the relative pose & the 'baseline score'
  C_Progress_display my_progress_bar( matches_provider_->pairWise_matches_.size(),
    std::cout,
    "Automatic selection of an initial pair:\n" );
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (const std::pair< Pair, IndMatches > & match_pair : matches_provider_->pairWise_matches_)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
#ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
#endif
      ++my_progress_bar;

      const Pair current_pair = match_pair.first;

      const size_t I = min(current_pair.first, current_pair.second);
      const size_t J = max(current_pair.first, current_pair.second);
      if (valid_views.count(I) && valid_views.count(J))
      {
        const View * view_I = sfm_data_.GetViews().at(I).get();
        const Intrinsics::const_iterator iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic);
        const View * view_J = sfm_data_.GetViews().at(J).get();
        const Intrinsics::const_iterator iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);

        const Pinhole_Intrinsic * cam_I = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic_I->second.get());
        const Pinhole_Intrinsic * cam_J = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic_J->second.get());
        if (cam_I != NULL && cam_J != NULL)
        {
          openMVG::tracks::STLMAPTracks map_tracksCommon;
          const std::set<size_t> set_imageIndex= {I, J};
          tracks::TracksUtilsMap::GetTracksInImages(set_imageIndex, map_tracks_, map_tracksCommon);

          // Copy points correspondences to arrays for relative pose estimation
          const size_t n = map_tracksCommon.size();
          Mat xI(2,n), xJ(2,n);
          size_t cptIndex = 0;
          for (openMVG::tracks::STLMAPTracks::const_iterator
            iterT = map_tracksCommon.begin(); iterT != map_tracksCommon.end();
            ++iterT, ++cptIndex)
          {
            tracks::submapTrack::const_iterator iter = iterT->second.begin();
            const size_t i = iter->second;
            const size_t j = (++iter)->second;

            Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
            xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
            feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
            xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
          }

          // Robust estimation of the relative pose
          RelativePose_Info relativePose_info;
          relativePose_info.initial_residual_tolerance = Square(4.0);

          if (robustRelativePose(
            cam_I->K(), cam_J->K(),
            xI, xJ, relativePose_info,
            std::make_pair(cam_I->w(), cam_I->h()), std::make_pair(cam_J->w(), cam_J->h()),
            256) && relativePose_info.vec_inliers.size() > iMin_inliers_count)
          {
            // Triangulate inliers & compute angle between bearing vectors
            std::vector<float> vec_angles;
            vec_angles.reserve(relativePose_info.vec_inliers.size());
            const Pose3 pose_I = Pose3(Mat3::Identity(), Vec3::Zero());
            const Pose3 pose_J = relativePose_info.relativePose;
            const Mat34 PI = cam_I->get_projective_equivalent(pose_I);
            const Mat34 PJ = cam_J->get_projective_equivalent(pose_J);
            for (const size_t inlier_idx : relativePose_info.vec_inliers)
            {
              Vec3 X;
              TriangulateDLT(PI, xI.col(inlier_idx), PJ, xJ.col(inlier_idx), &X);

              openMVG::tracks::STLMAPTracks::const_iterator iterT = map_tracksCommon.begin();
              std::advance(iterT, inlier_idx);
              tracks::submapTrack::const_iterator iter = iterT->second.begin();
              const Vec2 featI = features_provider_->feats_per_view[I][iter->second].coords().cast<double>();
              const Vec2 featJ = features_provider_->feats_per_view[J][(++iter)->second].coords().cast<double>();
              vec_angles.push_back(AngleBetweenRay(pose_I, cam_I, pose_J, cam_J, featI, featJ));
            }
            // Compute the median triangulation angle
            const unsigned median_index = vec_angles.size() / 2;
            std::nth_element(
              vec_angles.begin(),
              vec_angles.begin() + median_index,
              vec_angles.end());
            const float scoring_angle = vec_angles[median_index];
            // Store the pair iff the pair is in the asked angle range [fRequired_min_angle;fLimit_max_angle]
            if (scoring_angle > fRequired_min_angle &&
                scoring_angle < fLimit_max_angle)
            {
  #ifdef OPENMVG_USE_OPENMP
              #pragma omp critical
  #endif
              scoring_per_pair.emplace_back(scoring_angle, current_pair);
            }
          }
        }
      }
    } // omp section
  }
  std::sort(scoring_per_pair.begin(), scoring_per_pair.end());
  // Since scoring is ordered in increasing order, reverse the order
  std::reverse(scoring_per_pair.begin(), scoring_per_pair.end());
  if (!scoring_per_pair.empty())
  {
    initial_pair = scoring_per_pair.begin()->second;
    return true;
  }
  return false;
}

/// Compute the initial 3D seed (First camera t=0; R=Id, second estimated by 5 point algorithm)
bool IncrementalSfMReconstructionEngine::MakeInitialPair3D(const Pair & current_pair)
{
  // Compute robust Essential matrix for ImageId [I,J]
  // use min max to have I < J
  const size_t I = min(current_pair.first, current_pair.second);
  const size_t J = max(current_pair.first, current_pair.second);

  // a. Assert we have valid pinhole cameras
  const View * view_I = sfm_data_.GetViews().at(I).get();
  const Intrinsics::const_iterator iterIntrinsic_I = sfm_data_.GetIntrinsics().find(view_I->id_intrinsic);
  const View * view_J = sfm_data_.GetViews().at(J).get();
  const Intrinsics::const_iterator iterIntrinsic_J = sfm_data_.GetIntrinsics().find(view_J->id_intrinsic);

  if (iterIntrinsic_I == sfm_data_.GetIntrinsics().end() ||
      iterIntrinsic_J == sfm_data_.GetIntrinsics().end() )
  {
    return false;
  }

  const Pinhole_Intrinsic * cam_I = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic_I->second.get());
  const Pinhole_Intrinsic * cam_J = dynamic_cast<const Pinhole_Intrinsic*>(iterIntrinsic_J->second.get());
  if (cam_I == NULL || cam_J == NULL)
  {
    return false;
  }

  // b. Get common features between the two view
  // use the track to have a more dense match correspondence set
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  const std::set<size_t> set_imageIndex= {I, J};
  tracks::TracksUtilsMap::GetTracksInImages(set_imageIndex, map_tracks_, map_tracksCommon);

  //-- Copy point to arrays
  const size_t n = map_tracksCommon.size();
  Mat xI(2,n), xJ(2,n);
  size_t cptIndex = 0;
  for (openMVG::tracks::STLMAPTracks::const_iterator
    iterT = map_tracksCommon.begin(); iterT != map_tracksCommon.end();
    ++iterT, ++cptIndex)
  {
    tracks::submapTrack::const_iterator iter = iterT->second.begin();
    const size_t i = iter->second;
    const size_t j = (++iter)->second;

    Vec2 feat = features_provider_->feats_per_view[I][i].coords().cast<double>();
    xI.col(cptIndex) = cam_I->get_ud_pixel(feat);
    feat = features_provider_->feats_per_view[J][j].coords().cast<double>();
    xJ.col(cptIndex) = cam_J->get_ud_pixel(feat);
  }

  // c. Robust estimation of the relative pose
  RelativePose_Info relativePose_info;

  const std::pair<size_t, size_t> imageSize_I(cam_I->w(), cam_I->h());
  const std::pair<size_t, size_t> imageSize_J(cam_J->w(), cam_J->h());

  if (!robustRelativePose(
    cam_I->K(), cam_J->K(), xI, xJ, relativePose_info, imageSize_I, imageSize_J, 4096))
  {
    std::cerr << " /!\\ Robust estimation failed to compute E for this pair"
      << std::endl;
    return false;
  }
  std::cout << "A-Contrario initial pair residual: "
    << relativePose_info.found_residual_precision << std::endl;
  // Bound min precision at 1 pix.
  relativePose_info.found_residual_precision = std::max(relativePose_info.found_residual_precision, 1.0);

  map_ACThreshold_.insert(std::make_pair(I, relativePose_info.found_residual_precision));
  map_ACThreshold_.insert(std::make_pair(J, relativePose_info.found_residual_precision));
  set_remaining_view_id_.erase(view_I->id_view);
  set_remaining_view_id_.erase(view_J->id_view);

  // Log data for incremental SfM
  increment_camera_.first.emplace(view_I->id_view);
  increment_camera_.first.emplace(view_J->id_view);
  // Save parent of the added camera to map
  slam_pp_data.parent_cam_id[view_J->id_view] = std::numeric_limits<IndexT>::max();
  // Add parent information to the graph (second camera is the parent of first)
  slam_pp_data.addCamWithParentToGraph(view_J->id_view,std::numeric_limits<IndexT>::max());
  
  // Save parent of the added camera to map
  slam_pp_data.parent_cam_id[view_I->id_view] = view_J->id_view;
  slam_pp_data.addCamWithParentToGraph(view_I->id_view,view_J->id_view);

  // Save view info if we want to do order-style reconstruction
  if (bOrderedProcessing_)
  {
    lastUsedViewId_ = I;
  }

  const bool bRefine_using_BA = performInitialTwoViewBA_;

  if (bRefine_using_BA)
  {
    // Refine the defined scene
    SfM_Data tiny_scene;
    tiny_scene.views.insert(*sfm_data_.GetViews().find(view_I->id_view));
    tiny_scene.views.insert(*sfm_data_.GetViews().find(view_J->id_view));
    tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_I->id_intrinsic));
    tiny_scene.intrinsics.insert(*sfm_data_.GetIntrinsics().find(view_J->id_intrinsic));

    // Init poses
    const Pose3 & Pose_I = tiny_scene.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
    const Pose3 & Pose_J = tiny_scene.poses[view_J->id_pose] = relativePose_info.relativePose;

    // Init structure
    const Mat34 P1 = cam_I->get_projective_equivalent(Pose_I);
    const Mat34 P2 = cam_J->get_projective_equivalent(Pose_J);
    Landmarks & landmarks = tiny_scene.structure;

    for (openMVG::tracks::STLMAPTracks::const_iterator
      iterT = map_tracksCommon.begin();
      iterT != map_tracksCommon.end();
      ++iterT)
    {
      // Get corresponding points
      tracks::submapTrack::const_iterator iter = iterT->second.begin();
      const size_t i = iter->second;
      const size_t j = (++iter)->second;

      const Vec2 x1_ = features_provider_->feats_per_view[I][i].coords().cast<double>();
      const Vec2 x2_ = features_provider_->feats_per_view[J][j].coords().cast<double>();

      Vec3 X;
      TriangulateDLT(P1, x1_, P2, x2_, &X);
      Observations obs;
      obs[view_I->id_view] = Observation(x1_, i);
      obs[view_J->id_view] = Observation(x2_, j);
      landmarks[iterT->first].obs = std::move(obs);
      landmarks[iterT->first].X = X;
    }
    Save(tiny_scene, stlplus::create_filespec(sOut_directory_, "initialPair.ply"), ESfM_Data(ALL));

 
    // - refine only Structure and Rotations & translations (keep intrinsic constant)
    Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
    Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
    if (!bundle_adjustment_obj.Adjust(tiny_scene,
        Optimize_Options
        (
          Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
          Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
          Structure_Parameter_Type::ADJUST_ALL) // Adjust structure
        )
      )
    {
      return false;
    }

    // Save computed data
    const Pose3 pose_I = sfm_data_.poses[view_I->id_pose] = tiny_scene.poses[view_I->id_pose];
    const Pose3 pose_J = sfm_data_.poses[view_J->id_pose] = tiny_scene.poses[view_J->id_pose];

    // List inliers and save them
    for (Landmarks::const_iterator iter = tiny_scene.GetLandmarks().begin();
      iter != tiny_scene.GetLandmarks().end(); ++iter)
    {
      const IndexT trackId = iter->first;
      const Landmark & landmark = iter->second;
      const Observations & obs = landmark.obs;
      Observations::const_iterator iterObs_xI = obs.find(view_I->id_view);
      Observations::const_iterator iterObs_xJ = obs.find(view_J->id_view);

      const Observation & ob_xI = iterObs_xI->second;
      const Observation & ob_xJ = iterObs_xJ->second;

      const double angle = AngleBetweenRay(
        pose_I, cam_I, pose_J, cam_J, ob_xI.x, ob_xJ.x);
      const Vec2 residual_I = cam_I->residual(pose_I, landmark.X, ob_xI.x);
      const Vec2 residual_J = cam_J->residual(pose_J, landmark.X, ob_xJ.x);
      if ( angle > 2.0 &&
           pose_I.depth(landmark.X) > 0 &&
           pose_J.depth(landmark.X) > 0 &&
           residual_I.norm() < relativePose_info.found_residual_precision &&
           residual_J.norm() < relativePose_info.found_residual_precision)
      {
        sfm_data_.structure[trackId] = landmarks[trackId];
        
        // Log data for incremental SfM
        increment_structure_.first.emplace(trackId);
        (increment_observation_.first)[view_I->id_view].emplace(trackId);
        (increment_observation_.first)[view_J->id_view].emplace(trackId);
        slam_pp_data.owner_track_cam_id[trackId] = view_J->id_view;

      }
    }
  }
  else
  {
    // With no BA -> just add the points as they are triangulated
    // Init poses
    const Pose3 & pose_I = sfm_data_.poses[view_I->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
    const Pose3 & pose_J = sfm_data_.poses[view_J->id_pose] = relativePose_info.relativePose;

    // Init structure
    const Mat34 P1 = cam_I->get_projective_equivalent(pose_I);
    const Mat34 P2 = cam_J->get_projective_equivalent(pose_J);
    
    for (openMVG::tracks::STLMAPTracks::const_iterator
      iterT = map_tracksCommon.begin();
      iterT != map_tracksCommon.end();
      ++iterT)
    {
      // Get corresponding points
      tracks::submapTrack::const_iterator iter = iterT->second.begin();
      const size_t i = iter->second;
      const size_t j = (++iter)->second;

      const Vec2 x1_ = features_provider_->feats_per_view[I][i].coords().cast<double>();
      const Vec2 x2_ = features_provider_->feats_per_view[J][j].coords().cast<double>();

      Vec3 X;
      TriangulateDLT(P1, x1_, P2, x2_, &X);
      
      const Observation ob_xI = Observation(x1_, i);
      const Observation ob_xJ = Observation(x2_, j);

      const double angle = AngleBetweenRay(
        pose_I, cam_I, pose_J, cam_J, ob_xI.x, ob_xJ.x);
      const Vec2 residual_I = cam_I->residual(pose_I, X, ob_xI.x);
      const Vec2 residual_J = cam_J->residual(pose_J, X, ob_xJ.x);
      if ( angle > 2.0 &&
           pose_I.depth(X) > 0 &&
           pose_J.depth(X) > 0 &&
           residual_I.norm() < relativePose_info.found_residual_precision &&
           residual_J.norm() < relativePose_info.found_residual_precision)
      {
        Observations obs;
        obs[I] = Observation(x1_, i);
        obs[J] = Observation(x2_, j);
        sfm_data_.structure[iterT->first].obs = std::move(obs);
        sfm_data_.structure[iterT->first].X = X;

        // Log data for incremental SfM
        increment_structure_.first.emplace(iterT->first);
        (increment_observation_.first)[I].emplace(iterT->first);
        (increment_observation_.first)[J].emplace(iterT->first);
        slam_pp_data.owner_track_cam_id[iterT->first] =J;

      }
    }

    Save(sfm_data_, stlplus::create_filespec(sOut_directory_, "initialPair.ply"), ESfM_Data(ALL));
  }

  // Save outlier residual information
  Histogram<double> histoResiduals;
  std::cout << std::endl
    << "=========================\n"
    << " MSE Residual InitialPair Inlier: " << ComputeResidualsHistogram(&histoResiduals) << "\n"
    << "=========================" << std::endl;

  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    html_doc_stream_->pushInfo(htmlMarkup("h1","Essential Matrix."));
    ostringstream os;
    os << std::endl
      << "-------------------------------" << "<br>"
      << "-- Robust Essential matrix: <"  << I << "," <<J << "> images: "
      << view_I->s_Img_path << ","
      << view_J->s_Img_path << "<br>"
      << "-- Threshold: " << relativePose_info.found_residual_precision << "<br>"
      << "-- Resection status: " << "OK" << "<br>"
      << "-- Nb points used for robust Essential matrix estimation: "
      << xI.cols() << "<br>"
      << "-- Nb points validated by robust estimation: "
      << sfm_data_.structure.size() << "<br>"
      << "-- % points validated: "
      << sfm_data_.structure.size()/static_cast<float>(xI.cols())
      << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());

    html_doc_stream_->pushInfo(htmlMarkup("h2",
      "Residual of the robust estimation (Initial triangulation). Thresholded at: "
      + toString(relativePose_info.found_residual_precision)));

    html_doc_stream_->pushInfo(htmlMarkup("h2","Histogram of residuals"));

    std::vector<double> xBin = histoResiduals.GetXbinsValue();
    std::pair< std::pair<double,double>, std::pair<double,double> > range =
      autoJSXGraphViewport<double>(xBin, histoResiduals.GetHist());

    htmlDocument::JSXGraphWrapper jsxGraph;
    jsxGraph.init("InitialPairTriangulationKeptInfo",600,300);
    jsxGraph.addXYChart(xBin, histoResiduals.GetHist(), "line,point");
    jsxGraph.addLine(relativePose_info.found_residual_precision, 0,
      relativePose_info.found_residual_precision, histoResiduals.GetHist().front());
    jsxGraph.UnsuspendUpdate();
    jsxGraph.setViewport(range);
    jsxGraph.close();
    html_doc_stream_->pushInfo(jsxGraph.toStr());

    html_doc_stream_->pushInfo("<hr>");

    ofstream htmlFileStream( string(stlplus::folder_append_separator(sOut_directory_) +
      "Reconstruction_Report.html").c_str());
    htmlFileStream << html_doc_stream_->getDoc();
  }

  return !sfm_data_.structure.empty();
}

double IncrementalSfMReconstructionEngine::ComputeResidualsHistogram(Histogram<double> * histo)
{
  // Collect residuals for each observation
  std::vector<float> vec_residuals;
  vec_residuals.reserve(sfm_data_.structure.size());
  for(Landmarks::const_iterator iterTracks = sfm_data_.GetLandmarks().begin();
      iterTracks != sfm_data_.GetLandmarks().end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    for(Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const View * view = sfm_data_.GetViews().find(itObs->first)->second.get();
      const Pose3 pose = sfm_data_.GetPoseOrDie(view);
      const std::shared_ptr<IntrinsicBase> intrinsic = sfm_data_.GetIntrinsics().find(view->id_intrinsic)->second;
      const Vec2 residual = intrinsic->residual(pose, iterTracks->second.X, itObs->second.x);
      vec_residuals.push_back( fabs(residual(0)) );
      vec_residuals.push_back( fabs(residual(1)) );
    }
  }
  // Display statistics
  if (vec_residuals.size() > 1)
  {
    float dMin, dMax, dMean, dMedian;
    minMaxMeanMedian<float>(vec_residuals.begin(), vec_residuals.end(),
                            dMin, dMax, dMean, dMedian);
    if (histo)  {
      *histo = Histogram<double>(dMin, dMax, 10);
      histo->Add(vec_residuals.begin(), vec_residuals.end());
    }

    std::cout << std::endl << std::endl;
    std::cout << std::endl
      << "SequentialSfMReconstructionEngine::ComputeResidualsMSE." << "\n"
      << "\t-- #Tracks:\t" << sfm_data_.GetLandmarks().size() << std::endl
      << "\t-- Residual min:\t" << dMin << std::endl
      << "\t-- Residual median:\t" << dMedian << std::endl
      << "\t-- Residual max:\t "  << dMax << std::endl
      << "\t-- Residual mean:\t " << dMean << std::endl;

    return dMean;
  }
  return -1.0;
}

/// Functor to sort a vector of pair given the pair's second value
template<class T1, class T2, class Pred = std::less<T2> >
struct sort_pair_second {
  bool operator()(const std::pair<T1,T2>&left,
                    const std::pair<T1,T2>&right)
  {
    Pred p;
    return p(left.second, right.second);
  }
};

/**
 * @brief Estimate images on which we can compute the resectioning safely.
 *
 * @param[out] vec_possible_indexes: list of indexes we can use for resectioning.
 * @return False if there is no possible resection.
 *
 * Sort the images by the number of features id shared with the reconstruction.
 * Select the image I that share the most of correspondences.
 * Then keep all the images that have at least:
 *  0.75 * #correspondences(I) common correspondences to the reconstruction.
 */
bool IncrementalSfMReconstructionEngine::FindImagesWithPossibleResection(
  std::vector<size_t> & vec_possible_indexes)
{
  // Threshold used to select the best images
  static const float dThresholdGroup = dBAThresholdGroup_;

  vec_possible_indexes.clear();

  if (set_remaining_view_id_.empty() || sfm_data_.GetLandmarks().empty())
    return false;


  // Collect tracksIds
  std::set<size_t> reconstructed_trackId;
  Pair_Vec vec_putative; // ImageId, NbPutativeCommonPoint

  std::transform(sfm_data_.GetLandmarks().begin(), sfm_data_.GetLandmarks().end(),
    std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
    stl::RetrieveKey());


  if (bOrderedProcessing_)
  {
    std::set<size_t> viewsInWindow;

    while (viewsInWindow.empty() && lastUsedViewId_ < sfm_data_.views.size())
    {
      std::cout<<"\nFinding next candidate views (window "<<orderWindowSize_<<") from view "<<lastUsedViewId_<<"\n\n";
      for (std::set<size_t>::iterator it_set = set_remaining_view_id_.begin(); it_set != set_remaining_view_id_.end(); ++it_set)
      {
        const size_t v_id = *it_set;
        const size_t min_id = std::max<size_t>(0, lastUsedViewId_ - orderWindowSize_);
        const size_t max_id = std::min<size_t>(sfm_data_.views.size()-1, lastUsedViewId_ + orderWindowSize_);
        if (v_id >= min_id && v_id <= max_id)
        {
          viewsInWindow.emplace(v_id);
          std::cout<<" "<<v_id;
        }
      }
      if (viewsInWindow.empty())
      {
        lastUsedViewId_++;
      }
    }
    // If possible views were found check through the whole set
    if (viewsInWindow.empty())
    {
      if (bTryAllViews_)
      {
        std::cout<<"Force from rest of the views\n\n";
        // We try the standard method
      }
      else
      {
        // if we only want in processing in order we stop
        set_remaining_view_id_.clear();
        return false;
      }
    }
    else
    {
      // From the candidates in the window order by the number of matches
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
      for (std::set<size_t>::const_iterator iter = viewsInWindow.begin();
            iter != viewsInWindow.end(); ++iter)
      {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
        {
          const size_t viewId = *iter;

          // Compute 2D - 3D possible content
          openMVG::tracks::STLMAPTracks map_tracksCommon;
          const std::set<size_t> set_viewId = {viewId};
          tracks::TracksUtilsMap::GetTracksInImages(set_viewId, map_tracks_, map_tracksCommon);

          if (!map_tracksCommon.empty())
          {
            std::set<size_t> set_tracksIds;
            tracks::TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

            // Count the common possible putative point
            //  with the already 3D reconstructed trackId
            std::vector<size_t> vec_trackIdForResection;
            std::set_intersection(set_tracksIds.begin(), set_tracksIds.end(),
              reconstructed_trackId.begin(),
              reconstructed_trackId.end(),
              std::back_inserter(vec_trackIdForResection));

    #ifdef OPENMVG_USE_OPENMP
            #pragma omp critical
    #endif
            {
              vec_putative.push_back( make_pair(viewId, vec_trackIdForResection.size()));
            }
          }
        }
      }
      // Sort by the number of matches to the 3D scene.
      std::sort(vec_putative.begin(), vec_putative.end(), sort_pair_second<size_t, size_t, std::greater<size_t> >());

      std::cout<<"PUTATIVE: "<<vec_putative.size()<<"\n";

      for (size_t i = 0; i < vec_putative.size(); ++i)
      {
        vec_possible_indexes.push_back(vec_putative[i].first);
        std::cout<<"AA: "<<vec_possible_indexes[i]<<" :: "<<vec_putative[i].first<<"\n";
      }
      std::cout<<"POSSIBLE: "<<vec_possible_indexes.size()<<"\n";
      return true;
    }
  }

/*
  // Go through remaining images and add the ones fitting the window
  if (bOrderedProcessing_)
  {
    while (vec_possible_indexes.empty() && lastUsedViewId_ < sfm_data_.views.size())
    {
      std::cout<<"\nFinding next candidate views (window "<<orderWindowSize_<<") from view "<<lastUsedViewId_<<":\n";
      for (std::set<size_t>::iterator it_set = set_remaining_view_id_.begin(); it_set != set_remaining_view_id_.end(); ++it_set)
      {
        const size_t v_id = *it_set;
        const size_t min_id = std::max<size_t>(0, lastUsedViewId_ - orderWindowSize_);
        const size_t max_id = std::min<size_t>(sfm_data_.views.size()-1, lastUsedViewId_ + orderWindowSize_);
        if (v_id >= min_id && v_id <= max_id)
        {
          vec_possible_indexes.emplace_back(v_id);
          std::cout<<" "<<v_id;
        }
      }
      if (vec_possible_indexes.empty())
      {
        lastUsedViewId_++;
      }
    }
    // If possible views were found check through the whole set
    if (vec_possible_indexes.empty())
    {
      if (bTryAllViews_)
      {
        std::cout<<"Force from rest of the views\n\n";  

      }
      else
      {
        set_remaining_view_id_.clear();
        return false;
      }
    }
    else
    {
      return true;
    }
  }
*/


  // Find best view from all un-used views

  // Collect tracksIds
/* std::set<size_t> reconstructed_trackId;
  std::transform(sfm_data_.GetLandmarks().begin(), sfm_data_.GetLandmarks().end(),
    std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
    stl::RetrieveKey());

  Pair_Vec vec_putative; // ImageId, NbPutativeCommonPoint*/
#ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel
#endif
  for (std::set<size_t>::const_iterator iter = set_remaining_view_id_.begin();
        iter != set_remaining_view_id_.end(); ++iter)
  {
#ifdef OPENMVG_USE_OPENMP
  #pragma omp single nowait
#endif
    {
      const size_t viewId = *iter;

      // Compute 2D - 3D possible content
      openMVG::tracks::STLMAPTracks map_tracksCommon;
      const std::set<size_t> set_viewId = {viewId};
      tracks::TracksUtilsMap::GetTracksInImages(set_viewId, map_tracks_, map_tracksCommon);

      if (!map_tracksCommon.empty())
      {
        std::set<size_t> set_tracksIds;
        tracks::TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

        // Count the common possible putative point
        //  with the already 3D reconstructed trackId
        std::vector<size_t> vec_trackIdForResection;
        std::set_intersection(set_tracksIds.begin(), set_tracksIds.end(),
          reconstructed_trackId.begin(),
          reconstructed_trackId.end(),
          std::back_inserter(vec_trackIdForResection));

#ifdef OPENMVG_USE_OPENMP
        #pragma omp critical
#endif
        {
          vec_putative.push_back( make_pair(viewId, vec_trackIdForResection.size()));
        }
      }
    }
  }

  // Sort by the number of matches to the 3D scene.
  std::sort(vec_putative.begin(), vec_putative.end(), sort_pair_second<size_t, size_t, std::greater<size_t> >());

  // If the list is empty or if the list contains images with no correspdences
  // -> (no resection will be possible)
  if (vec_putative.empty() || vec_putative[0].second == 0)
  {
    // All remaining images cannot be used for pose estimation
    set_remaining_view_id_.clear();
    return false;
  }

  // Add the image view index that share the most of 2D-3D correspondences
  vec_possible_indexes.push_back(vec_putative[0].first);

  int added_points = 0;

  // Then, add all the image view indexes that have at least N% of the number of the matches of the best image.
  const IndexT M = vec_putative[0].second; // Number of 2D-3D correspondences
  const size_t threshold = static_cast<size_t>(dThresholdGroup * M);
  for (size_t i = 1; i < vec_putative.size() &&
    vec_putative[i].second > threshold; ++i)
  {
    vec_possible_indexes.push_back(vec_putative[i].first);
    added_points+=vec_putative[i].second;
  }
  
  return true;
}






/**
 * @brief Add one image to the 3D reconstruction. To the resectioning of
 * the camera and triangulate all the new possible tracks.
 * @param[in] viewIndex: image index to add to the reconstruction.
 *
 * A. Compute 2D/3D matches
 * B. Look if intrinsic data is known or not
 * C. Do the resectioning: compute the camera pose.
 * D. Refine the pose of the found camera
 * E. Update the global scene with the new camera
 * F. Update the observations into the global scene structure
 * G. Triangulate new possible 2D tracks
 */
bool IncrementalSfMReconstructionEngine::Resection(const size_t viewIndex)
{
  using namespace tracks;

  // A. Compute 2D/3D matches
  // A1. list tracks ids used by the view
  openMVG::tracks::STLMAPTracks map_tracksCommon;
  const std::set<size_t> set_viewIndex = {viewIndex};
  TracksUtilsMap::GetTracksInImages(set_viewIndex, map_tracks_, map_tracksCommon);
  std::set<size_t> set_tracksIds;
  TracksUtilsMap::GetTracksIdVector(map_tracksCommon, &set_tracksIds);

  // A2. intersects the track list with the reconstructed
  std::set<size_t> reconstructed_trackId;
  std::transform(sfm_data_.GetLandmarks().begin(), sfm_data_.GetLandmarks().end(),
    std::inserter(reconstructed_trackId, reconstructed_trackId.begin()),
    stl::RetrieveKey());

  // Get the ids of the already reconstructed tracks
  std::set<size_t> set_trackIdForResection;
  std::set_intersection(set_tracksIds.begin(), set_tracksIds.end(),
    reconstructed_trackId.begin(),
    reconstructed_trackId.end(),
    std::inserter(set_trackIdForResection, set_trackIdForResection.begin()));

  if (set_trackIdForResection.empty())
  {
    // No match. The image has no connection with already reconstructed points.
    std::cout << std::endl
      << "-------------------------------" << "\n"
      << "-- Resection of camera index: " << viewIndex << "\n"
      << "-- Resection status: " << "FAILED" << "\n"
      << "-------------------------------" << std::endl;
    return false;
  }

  // Get back featId associated to a tracksID already reconstructed.
  // These 2D/3D associations will be used for the resection.
  std::vector<size_t> vec_featIdForResection;
  TracksUtilsMap::GetFeatIndexPerViewAndTrackId(map_tracksCommon,
    set_trackIdForResection,
    viewIndex,
    &vec_featIdForResection);

  // Localize the image inside the SfM reconstruction
  Image_Localizer_Match_Data resection_data;
  resection_data.pt2D.resize(2, set_trackIdForResection.size());
  resection_data.pt3D.resize(3, set_trackIdForResection.size());

  // B. Look if intrinsic data is known or not
  const View * view_I = sfm_data_.GetViews().at(viewIndex).get();
  std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic (nullptr);
  if (sfm_data_.GetIntrinsics().count(view_I->id_intrinsic))
  {
    optional_intrinsic = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic);
  }

  // Setup the track 2d observation for this new view
  Mat2X pt2D_original(2, set_trackIdForResection.size());
  std::set<size_t>::const_iterator iterTrackId = set_trackIdForResection.begin();
  std::vector<size_t>::const_iterator iterfeatId = vec_featIdForResection.begin();
  for (size_t cpt = 0; cpt < vec_featIdForResection.size(); ++cpt, ++iterTrackId, ++iterfeatId)
  {
    resection_data.pt3D.col(cpt) = sfm_data_.GetLandmarks().at(*iterTrackId).X;
    resection_data.pt2D.col(cpt) = pt2D_original.col(cpt) =
      features_provider_->feats_per_view.at(viewIndex)[*iterfeatId].coords().cast<double>();
    // Handle image distortion if intrinsic is known (to ease the resection)
    if (optional_intrinsic && optional_intrinsic->have_disto())
    {
      resection_data.pt2D.col(cpt) = optional_intrinsic->get_ud_pixel(resection_data.pt2D.col(cpt));
    }
  }

  // C. Do the resectioning: compute the camera pose
  std::cout << std::endl
    << "-------------------------------" << std::endl
    << "-- Robust Resection of view: " << viewIndex << std::endl;

  geometry::Pose3 pose;
  const bool bResection = sfm::SfM_Localizer::Localize
  (
    Pair(view_I->ui_width, view_I->ui_height),
    optional_intrinsic.get(),
    resection_data,
    pose
  );
  resection_data.pt2D = std::move(pt2D_original); // restore original image domain points

  if (!sLogging_file_.empty())
  {
    using namespace htmlDocument;
    ostringstream os;
    os << "Resection of Image index: <" << viewIndex << "> image: "
      << view_I->s_Img_path <<"<br> \n";
    html_doc_stream_->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << std::endl
      << "-------------------------------" << "<br>"
      << "-- Robust Resection of camera index: <" << viewIndex << "> image: "
      <<  view_I->s_Img_path <<"<br>"
      << "-- Threshold: " << resection_data.error_max << "<br>"
      << "-- Resection status: " << (bResection ? "OK" : "FAILED") << "<br>"
      << "-- Nb points used for Resection: " << vec_featIdForResection.size() << "<br>"
      << "-- Nb points validated by robust estimation: " << resection_data.vec_inliers.size() << "<br>"
      << "-- % points validated: "
      << resection_data.vec_inliers.size()/static_cast<float>(vec_featIdForResection.size()) << "<br>"
      << "-------------------------------" << "<br>";
    html_doc_stream_->pushInfo(os.str());
  }

  if (!bResection)
    return false;

  // D. Refine the pose of the found camera.
  // We use a local scene with only the 3D points and the new camera.
  {
    const bool b_new_intrinsic = (optional_intrinsic == nullptr);
    // A valid pose has been found (try to refine it):
    // If no valid intrinsic as input:
    //  init a new one from the projection matrix decomposition
    // Else use the existing one and consider it as constant.
    if (b_new_intrinsic)
    {
      // setup a default camera model from the found projection matrix
      Mat3 K, R;
      Vec3 t;
      KRt_From_P(resection_data.projection_matrix, &K, &R, &t);

      const double focal = (K(0,0) + K(1,1))/2.0;
      const Vec2 principal_point(K(0,2), K(1,2));

      // Create the new camera intrinsic group
      switch (cam_type_)
      {
        case PINHOLE_CAMERA:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_RADIAL1:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic_Radial_K1>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_RADIAL3:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic_Radial_K3>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_BROWN:
          optional_intrinsic =
            std::make_shared<Pinhole_Intrinsic_Brown_T2>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        case PINHOLE_CAMERA_FISHEYE:
            optional_intrinsic =
                std::make_shared<Pinhole_Intrinsic_Fisheye>
            (view_I->ui_width, view_I->ui_height, focal, principal_point(0), principal_point(1));
        break;
        default:
          std::cerr << "Try to create an unknown camera type." << std::endl;
          return false;
      }
    }
    if(performLocalPoseBA_)
    {
      const bool b_refine_pose = true;
      const bool b_refine_intrinsics = false;
      if(!sfm::SfM_Localizer::RefinePose(
          optional_intrinsic.get(), pose,
          resection_data, b_refine_pose, b_refine_intrinsics))
      {
        return false;
      }
    }

    // E. Update the global scene with:
    // - the new found camera pose
    sfm_data_.poses[view_I->id_pose] = pose;
    // - track the view's AContrario robust estimation found threshold
    map_ACThreshold_.insert(std::make_pair(viewIndex, resection_data.error_max));
    // - intrinsic parameters (if the view has no intrinsic group add a new one)
    if (b_new_intrinsic)
    {
      // Since the view have not yet an intrinsic group before, create a new one
      IndexT new_intrinsic_id = 0;
      if (!sfm_data_.GetIntrinsics().empty())
      {
        // Since some intrinsic Id already exists,
        //  we have to create a new unique identifier following the existing one
        std::set<IndexT> existing_intrinsicId;
          std::transform(sfm_data_.GetIntrinsics().begin(), sfm_data_.GetIntrinsics().end(),
          std::inserter(existing_intrinsicId, existing_intrinsicId.begin()),
          stl::RetrieveKey());
        new_intrinsic_id = (*existing_intrinsicId.rbegin())+1;
      }
      sfm_data_.views.at(viewIndex).get()->id_intrinsic = new_intrinsic_id;
      sfm_data_.intrinsics[new_intrinsic_id] = optional_intrinsic;
    }



    // Local cams in common BA
    
    if(performLocalPoseBA_)
    {
    // we have the tracks that are seen by the new camera
    // Check which already reconstructed views see it
      {
        std::cout<<"Performing Local BA\n";
        std::set<IndexT> affected_views;
        std::set<IndexT> valid_views = Get_Valid_Views(sfm_data_);


        for (auto trackId: set_trackIdForResection)
        { 
          Landmark & landmark = sfm_data_.structure[trackId];
          for (auto obs_i : landmark.obs)
          {
            if (valid_views.find(obs_i.first) != valid_views.end())
            {
              affected_views.emplace(obs_i.first);
              // Add view to tiny scene if its not yet added
            }
          }
        }
        std::cout<<"Affected views: "<<affected_views.size()<<"\n";
        std::cout<<"Affected landmarks: "<<set_trackIdForResection.size()<<"\n";
        // Go throught all the tracks seen by new view and add them to the scene

        //Save(sfm_data_, stlplus::create_filespec(sOut_directory_, "scene_before_local_all.ply"), ESfM_Data(ALL));

        Bundle_Adjustment_Ceres::BA_Ceres_options options(true, true);
        options.bUse_loss_function_ = true;
        options.linear_solver_type_ = ceres::DENSE_SCHUR;
        Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
        bundle_adjustment_obj.LocalAdjust(sfm_data_,
            Optimize_Options
            (
              Intrinsic_Parameter_Type::NONE, // Keep intrinsic constant
              Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
              Structure_Parameter_Type::ADJUST_ALL // Adjust structure
            ),
            affected_views,
            set_trackIdForResection);

        //Save(sfm_data_, stlplus::create_filespec(sOut_directory_, "scene_after_local_all.ply"), ESfM_Data(ALL));
   
      }
    }
    // Log data for incremental SfM 
    increment_camera_.first.emplace(viewIndex);
    
    // Find parent camera (camera that gives most inliers for pose estimation)
    {
      Hash_Map<IndexT,std::set<IndexT> >resectionTracksPerView;
      // Get all reconstructed views
      std::set<IndexT> valid_views = Get_Valid_Views(sfm_data_);
      // remove current view from the list
      valid_views.erase(valid_views.find(view_I->id_pose));

      // iterate only over inlier tracks
      for (auto it_inlier : resection_data.vec_inliers)
      {
        std::set<size_t>::iterator iterTrackId = set_trackIdForResection.begin();
        std::advance(iterTrackId,it_inlier);
        size_t trackId = *iterTrackId;  
        for (const auto & iterT : map_tracks_[trackId])
        {
          if (valid_views.find(iterT.first) != valid_views.end())
          {
            resectionTracksPerView[iterT.first].emplace(trackId);
          }
        }
      }
      size_t max_n_tracks_cam = 0;
      IndexT parent_cam_id;
      for (auto it_r : resectionTracksPerView)
      {
        if (it_r.second.size() > max_n_tracks_cam)
        {
          max_n_tracks_cam = it_r.second.size();
          parent_cam_id = it_r.first;
        }
      }

      std::cout << "Camera: "<<view_I->id_pose<<" Parent camera: "<< parent_cam_id << " with " << max_n_tracks_cam<< " inliers\n";
      // Save parent of the added camera to map
      slam_pp_data.parent_cam_id[view_I->id_view] = parent_cam_id;
      slam_pp_data.addCamWithParentToGraph(view_I->id_pose, parent_cam_id);

    }
  }

  // F. List tracks that share content with this view and add observations and new 3D track if required.
  //    - If the track already exists (look if the new view tracks observation are valid)
  //    - If the track does not exists, try robust triangulation & add the new valid view track observation
  {
    // Get information of new view
    const IndexT I = viewIndex;
    const View * view_I = sfm_data_.GetViews().at(I).get();
    const IntrinsicBase * cam_I = sfm_data_.GetIntrinsics().at(view_I->id_intrinsic).get();
    const Pose3 pose_I = sfm_data_.GetPoseOrDie(view_I);

    // Vector of all already reconstructed views
    const std::set<IndexT> valid_views = Get_Valid_Views(sfm_data_);

    // Go through each track and look if we must add new view observations or new 3D points
    for (const std::pair< size_t, tracks::submapTrack >& trackIt : map_tracksCommon)
    {
      const size_t trackId = trackIt.first;
      const tracks::submapTrack & track = trackIt.second;

      // List the potential view observations of the track
      const tracks::submapTrack & allViews_of_track = map_tracks_[trackId];

      // List to save the new view observations that must be added to the track
      std::set<IndexT> new_track_observations_valid_views;

      // If the track was already reconstructed
      if (sfm_data_.structure.count(trackId) != 0)
      {
        // Since the 3D point was triangulated before we add the new the Inth view observation
        new_track_observations_valid_views.insert(I);
      }
      else
      {
        // Go through the views that observe this track & look if a successful triangulation can be done
        for (const std::pair< IndexT, IndexT >& trackViewIt : allViews_of_track)
        {
          const IndexT & J = trackViewIt.first;
          // If view is valid try triangulation
          if (J!=I && valid_views.count(J) != 0 )
          {
            // If successfuly triangulated add the observation from J view
            if (sfm_data_.structure.count(trackId) != 0)
            {
              new_track_observations_valid_views.insert(J);
            }
            else
            {
              const View * view_J = sfm_data_.GetViews().at(J).get();
              const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();
              const Pose3 pose_J = sfm_data_.GetPoseOrDie(view_J);
              const Vec2 xJ = features_provider_->feats_per_view.at(J)[allViews_of_track.at(J)].coords().cast<double>();

              // Position of the point in view I
              const Vec2 xI = features_provider_->feats_per_view.at(I)[track.at(I)].coords().cast<double>();

              // Try to triangulate a 3D point from J view
              // A new 3D point must be added
              // Triangulate it
              const Vec2 xI_ud = cam_I->get_ud_pixel(xI);
              const Vec2 xJ_ud = cam_J->get_ud_pixel(xJ);
              const Mat34 P_I = cam_I->get_projective_equivalent(pose_I);
              const Mat34 P_J = cam_J->get_projective_equivalent(pose_J);
              Vec3 X = Vec3::Zero();
              TriangulateDLT(P_I, xI_ud, P_J, xJ_ud, &X);
              // Check triangulation result
              const double angle = AngleBetweenRay(pose_I, cam_I, pose_J, cam_J, xI, xJ);
              const Vec2 residual_I = cam_I->residual(pose_I, X, xI);
              const Vec2 residual_J = cam_J->residual(pose_J, X, xJ);
              if (
                  //  - Check angle (small angle leads to imprecise triangulation)
                  angle > 2.0 &&
                  //  - Check positive depth
                  pose_I.depth(X) > 0 &&
                  pose_J.depth(X) > 0 &&
                  //  - Check residual values (must be inferior to the found view's AContrario threshold)
                  residual_I.norm() < std::max(4.0, map_ACThreshold_.at(I)) &&
                  residual_J.norm() < std::max(4.0, map_ACThreshold_.at(J))
                 )
              {
                // Add a new track
                Landmark & landmark = sfm_data_.structure[trackId];
                landmark.X = X;
                new_track_observations_valid_views.insert(I);
                new_track_observations_valid_views.insert(J);

                // Log data for incremental SfM 
                increment_structure_.first.emplace(trackId);
                slam_pp_data.owner_track_cam_id[trackId] = viewIndex;
              } // 3D point is valid
              else
              {
                // We mark the view to add the observations once the point is triangulated
                new_track_observations_valid_views.insert(J);
              } // 3D point is invalid
            }
          }
        }// Go through all the views
      }// If new point

      // If successfuly triangulated, add the valid view observations
      if (sfm_data_.structure.count(trackId) != 0 &&
          !new_track_observations_valid_views.empty()
         )
      {
        Landmark & landmark = sfm_data_.structure[trackId];
        // Check if view feature point observations of the track are valid (residual, depth) or not
        for (const IndexT &J: new_track_observations_valid_views)
        {
          const View * view_J = sfm_data_.GetViews().at(J).get();
          const IntrinsicBase * cam_J = sfm_data_.GetIntrinsics().at(view_J->id_intrinsic).get();
          const Pose3 pose_J = sfm_data_.GetPoseOrDie(view_J);
          const Vec2 xJ = features_provider_->feats_per_view.at(J)[allViews_of_track.at(J)].coords().cast<double>();

          // If we dont do local outlier removal we automatically add observation to the system without checking anything 
          if (performLocalOutlierRemoval_)
          {
            const Vec2 residual = cam_J->residual(pose_J, landmark.X, xJ);
            if (pose_J.depth(landmark.X) > 0 &&
                residual.norm() < std::max(4.0, map_ACThreshold_.at(J))
               )
            {
              landmark.obs[J] = Observation(xJ, allViews_of_track.at(J));

              // Log data for incremental SfM 
              (increment_observation_.first)[view_J->id_view].emplace(trackId);
            }
          }
          else
          {
            landmark.obs[J] = Observation(xJ, allViews_of_track.at(J));

            // Log data for incremental SfM 
            (increment_observation_.first)[view_J->id_view].emplace(trackId);
          }
        }
      }
    }// All the tracks in the view
  }

  return true;
}

/// Bundle adjustment to refine Structure; Motion and Intrinsics
bool IncrementalSfMReconstructionEngine::BundleAdjustment()
{
  Bundle_Adjustment_Ceres::BA_Ceres_options options;
  if ( sfm_data_.GetPoses().size() > 100 &&
      (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
       ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
       ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
      )
  // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
  {
    options.preconditioner_type_ = ceres::JACOBI;
    options.linear_solver_type_ = ceres::SPARSE_SCHUR;
  }
  else
  {
    options.linear_solver_type_ = ceres::DENSE_SCHUR;
  }
  Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
  const Optimize_Options ba_refine_options
    ( ReconstructionEngine::intrinsic_refinement_options_,
      Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
      Structure_Parameter_Type::ADJUST_ALL // Adjust scene structure
    );
  return bundle_adjustment_obj.Adjust(sfm_data_, ba_refine_options);
}

/**
 * @brief Discard tracks with too large residual error
 *
 * Remove observation/tracks that have:
 *  - too large residual error
 *  - too small angular value
 *
 * @return True if more than 'count' outliers have been removed.
 */
bool IncrementalSfMReconstructionEngine::badTrackRejector(double dPrecision, size_t count)
{
  // Funcions are added here to be able to keep track of removed objects
  //const size_t nbOutliers_residualErr = RemoveOutliers_PixelResidualError(sfm_data_, dPrecision, 2);
  double dThresholdPixel = dPrecision;
  unsigned int minTrackLength = 2;

  IndexT nbOutliers_residualErr = 0;
  Landmarks::iterator iterTracks = sfm_data_.structure.begin();
  while (iterTracks != sfm_data_.structure.end())
  {
    Observations & obs = iterTracks->second.obs;
    Observations::iterator itObs = obs.begin();
    bool dStructure=false;
    while (!dStructure && itObs != obs.end())
    {
      const View * view = sfm_data_.views.at(itObs->first).get();
      const geometry::Pose3 pose = sfm_data_.GetPoseOrDie(view);
      const cameras::IntrinsicBase * intrinsic = sfm_data_.intrinsics.at(view->id_intrinsic).get();
      const Vec2 residual = intrinsic->residual(pose, iterTracks->second.X, itObs->second.x);
      if (residual.norm() > dThresholdPixel)
      {
        // Point is rejected
        // Current HACK - we dont delete observation if its owner
        // Check if we are removing the owners measurement --> asign a different owner
        if ( slam_pp_data.owner_track_cam_id[iterTracks->first] != itObs->first)
        {
          // Mark that the observation was removed
          observation_last_removed[itObs->first][iterTracks->first] = current_recon_iteration;
          std::cout<<"Remove point as the core is wrong\n";
          obs.clear();
          dStructure = true;
          
          /*++nbOutliers_residualErr;
          itObs = obs.erase(itObs);*/
        }
      else
        ++itObs;
      }
      else
        ++itObs;
    }
    if (obs.empty() || obs.size() < minTrackLength)
    {
      // Mark that the structure was removed
      structure_last_removed[iterTracks->first] = current_recon_iteration;
      std::cout<<"Remove structure\n";
      iterTracks = sfm_data_.structure.erase(iterTracks);

    }
    else
      ++iterTracks;
  }

  // const size_t nbOutliers_angleErr = RemoveOutliers_AngleError(sfm_data_, 2.0);
  double dMinAcceptedAngle = 2.0;

  IndexT nbOutliers_angleErr = 0;
  iterTracks = sfm_data_.structure.begin();
  while (iterTracks != sfm_data_.structure.end())
  {
    Observations & obs = iterTracks->second.obs;
    double max_angle = 0.0;
    for (Observations::const_iterator itObs1 = obs.begin();
      itObs1 != obs.end(); ++itObs1)
    {
      const View * view1 = sfm_data_.views.at(itObs1->first).get();
      const geometry::Pose3 pose1 = sfm_data_.GetPoseOrDie(view1);
      const cameras::IntrinsicBase * intrinsic1 = sfm_data_.intrinsics.at(view1->id_intrinsic).get();

      Observations::const_iterator itObs2 = itObs1;
      ++itObs2;
      for (; itObs2 != obs.end(); ++itObs2)
      {
        const View * view2 = sfm_data_.views.at(itObs2->first).get();
        const geometry::Pose3 pose2 = sfm_data_.GetPoseOrDie(view2);
        const cameras::IntrinsicBase * intrinsic2 = sfm_data_.intrinsics.at(view2->id_intrinsic).get();

        const double angle = AngleBetweenRay(
          pose1, intrinsic1, pose2, intrinsic2,
          itObs1->second.x, itObs2->second.x);
        max_angle = std::max(angle, max_angle);
      }
    }
    if (max_angle < dMinAcceptedAngle)
    {
      // Mark that the structure was removed
      structure_last_removed[iterTracks->first] = current_recon_iteration;
  std::cout<<"Remove structure\n";
      iterTracks = sfm_data_.structure.erase(iterTracks);
      ++nbOutliers_angleErr;

    }
    else
      ++iterTracks;
  }

  return (nbOutliers_residualErr + nbOutliers_angleErr) > count;
}

/// Check if incremental history is consistent with the data we logged
bool IncrementalSfMReconstructionEngine::checkIncrementalConsistency()
{
  std::cout<<"START CONSISTENCY CHECK\n";
  std::vector<std::pair<std::set<IndexT>, std::set<IndexT> > > h_cam = history_i_camera_;
  std::vector<std::pair<std::set<IndexT>, std::set<IndexT> > > h_struct = history_i_structure_;
  std::vector<std::pair< Hash_Map<IndexT, std::set<IndexT> >, Hash_Map<IndexT, std::set<IndexT> > > > h_obs = history_i_observations_;


  size_t n_start_cams = 0;
  size_t n_start_struct = 0;
  size_t n_start_obs = 0;
  size_t n_err_cams = 0;
  size_t n_err_struct = 0;
  size_t n_err_obs = 0;
  size_t n_left_cams = 0;
  size_t n_left_struct = 0;
  size_t n_left_obs = 0;

  // Get initial data
  for (size_t c_i = 0; c_i<h_cam.size();c_i++)
  {
    n_start_cams += h_cam[c_i].first.size();
  }
  for (size_t c_i = 0; c_i<h_struct.size();c_i++)
  {
    n_start_struct += h_struct[c_i].first.size();
  }
  for (size_t c_i = 0; c_i<h_obs.size();c_i++)
  {
    Hash_Map<IndexT, std::set<IndexT> > & views_tracks = h_obs[c_i].first;
    for (Hash_Map<IndexT, std::set<IndexT> >::iterator it_view_tracks = views_tracks.begin(); it_view_tracks!= views_tracks.end(); ++it_view_tracks)
    {
      n_start_obs +=  it_view_tracks->second.size();
    }
  }


  // Check if all cameras are present
  const std::set<IndexT> valid_views = Get_Valid_Views(sfm_data_);
  for (std::set<IndexT>::iterator it_valid_v = valid_views.begin(); it_valid_v!=valid_views.end();++it_valid_v)
  {
    IndexT valid_view_id = *it_valid_v;
    bool found_element = false;
    IndexT c_i = 0;
    while (c_i < h_cam.size() && !found_element)
    {
      std::set<IndexT> & cams_in_step = h_cam[c_i].first;
      std::set<IndexT>::iterator it_find_cam = std::find(cams_in_step.begin(),cams_in_step.end(),valid_view_id);
      if (it_find_cam!=cams_in_step.end())
      {
        found_element = true;
        cams_in_step.erase(it_find_cam);
      }
      c_i++;
    }
    if (!found_element)
    {
      n_err_cams++;
      std::cerr<<"Camera "<<valid_view_id<<" Not found!!!\n";
    }
  }

  
  Landmarks & landmarks = sfm_data_.structure;
  for (Landmarks::iterator it_landmark = landmarks.begin(); it_landmark!=landmarks.end();++it_landmark)
  {
    IndexT landmark_id = it_landmark->first;
    Landmark & landmark = it_landmark->second;
    bool found_structure = false;
    IndexT s_i = 0;
    while (s_i < h_struct.size() && !found_structure)
    {
      std::set<IndexT> & structure_in_step = h_struct[s_i].first;
      std::set<IndexT>::iterator it_find_structure = std::find(structure_in_step.begin(),structure_in_step.end(),landmark_id);
      if (it_find_structure!=structure_in_step.end())
      {
        found_structure = true;
        structure_in_step.erase(it_find_structure);
        // Try to find observations
        for (Observations::iterator it_obs = landmark.obs.begin(); it_obs != landmark.obs.end(); ++it_obs)
        {
          IndexT view_id = it_obs->first;
          
          bool found_observation = false;
          IndexT o_i = 0;
          while (o_i < h_obs.size() && !found_observation)
          {
            Hash_Map<IndexT, std::set<IndexT> > & obs_in_step = h_obs[o_i].first;
            for (Hash_Map<IndexT, std::set<IndexT> >::iterator it_obs_step = obs_in_step.begin(); it_obs_step != obs_in_step.end(); ++it_obs_step)
            {
              if (it_obs_step->first == view_id)
              {
                std::set<IndexT> &tracks_id_per_view = it_obs_step->second;
                std::set<IndexT>::iterator it_track_view = std::find(tracks_id_per_view.begin(),tracks_id_per_view.end(),landmark_id);
                if (it_track_view != tracks_id_per_view.end())
                {
                  found_observation = true;
                  tracks_id_per_view.erase(it_track_view);
                  break;
                }
              }
            }
            o_i++;
          }
          if (!found_observation)
          {
            n_err_obs++;
            std::cerr<<"Observation "<<view_id<<" :: "<<landmark_id<<" Not found!!!\n";
          }
        }
      }
      s_i++;
    }
    if (!found_structure)
    {
      n_err_struct++;
      std::cerr<<"Landmark "<<landmark_id<<" Not found!!!\n";
    }
  }

  // Check if anything extra in incremental log
  for (size_t c_i = 0; c_i<h_cam.size();c_i++)
  {
    n_left_cams += h_cam[c_i].first.size();
  }
  for (size_t c_i = 0; c_i<h_struct.size();c_i++)
  {
    n_left_struct += h_struct[c_i].first.size();
  }
  for (size_t c_i = 0; c_i<h_obs.size();c_i++)
  {
    Hash_Map<IndexT, std::set<IndexT> > & views_tracks = h_obs[c_i].first;
    for (Hash_Map<IndexT, std::set<IndexT> >::iterator it_view_tracks = views_tracks.begin(); it_view_tracks!= views_tracks.end(); ++it_view_tracks)
    {
      n_left_obs +=  it_view_tracks->second.size();
    }
  }
  std::cout << "Cameras start / errors / left : " << n_start_cams << " / " << n_err_cams << " / " << n_left_cams << "\n";
  std::cout << "Landmarks start / errors / left : " << n_start_struct << " / " << n_err_struct << " / " << n_left_struct << "\n";
  std::cout << "Observation start / errors / left : " << n_start_obs << " / " << n_err_obs << " / " << n_left_obs << "\n";


  return true;
}

void IncrementalSfMReconstructionEngine::resetIncrementStep()
{
  // Save history (for debug)
  //if (performConsistencyCheck_)
  //{
    history_i_camera_.emplace_back(increment_camera_);
    history_i_structure_.emplace_back(increment_structure_);
    history_i_observations_.emplace_back(increment_observation_);
  //}
  // Clear for next step
  increment_camera_.first.clear();
  increment_structure_.first.clear();
  increment_observation_.first.clear();
}

/*
  void IncrementalSfMReconstructionEngine::ExportIncrementToGraphFile_SlamPP()
  {

    // Cameras
    for (std::set<IndexT>::iterator c_i = increment_camera_.first.begin(); c_i != increment_camera_.first.end(); ++c_i)
    {
      const IndexT camId_omvg = *c_i;
      const IndexT camId_slamPP = slam_pp_data.getNextFreeSlamPPId();

      // Save mapping between omvg and slamPP
      slam_pp_data.setCamId_SlamPP(camId_omvg,camId_slamPP);

      // Get information of the camera
      const View * view = sfm_data_.GetViews().at(camId_omvg).get();
      if ( !sfm_data_.IsPoseAndIntrinsicDefined( view ) )
      {
        continue;
      }

      const IntrinsicBase * cam = sfm_data_.GetIntrinsics().at(view->id_intrinsic).get();
      const Pinhole_Intrinsic * pinhole_cam = static_cast<const Pinhole_Intrinsic *>( cam );
      const Mat3 K = pinhole_cam->K();

      switch (slam_pp_data.iOutputVertexType)
      {
        case 0: // SE(3)
        {
          // Get pose
          Pose3 pose= sfm_data_.GetPoseOrDie(view);
          const Mat3 rotation = pose.rotation().transpose();
          const Vec3 center = pose.center();
          Eigen::Quaterniond q( rotation ) ;
          // Export to graph file
          slam_pp_data.slamPP_DatasetFile << "VERTEX_CAM" 
            << " " << camId_slamPP
            << " " << center[0]
            << " " << center[1]
            << " " << center[2]
            << " " << q.x()
            << " " << q.y()
            << " " << q.z()
            << " " << q.w()
            << " " << K(0,0)
            << " " << K(1,1)
            << " " << K(0,2)
            << " " << K(1,2)
            << " " << "0.0"
            << std::endl;
        }
        break;
        case 1: // Sim(3)
        {
          Pose3 pose = sfm_data_.GetPoseOrDie(view);
          const Mat3 rotation = pose.rotation().transpose();
          const Vec3 center = pose.center();
          const double scale = 1.0;
          double angleAxis[3];
          ceres::RotationMatrixToAngleAxis((const double*)rotation.data(), angleAxis);
          // Export to graph file
          slam_pp_data.slamPP_DatasetFile << "VERTEX_CAM:SIM3" 
            << " " << camId_slamPP
            << " " << center[0]
            << " " << center[1]
            << " " << center[2]
            << " " << angleAxis[0]
            << " " << angleAxis[1]
            << " " << angleAxis[2]
            << " " << scale
            << " " << K(0,0)
            << " " << K(1,1)
            << " " << K(0,2)
            << " " << K(1,2)
            << " " << "0.0"
            << std::endl;
        }
        break;
      }
    
    }

    // Structure
    Landmarks & landmarks = sfm_data_.structure;
    for (std::set<IndexT>::iterator s_i = increment_structure_.first.begin(); s_i != increment_structure_.first.end(); ++s_i)
    {
      const IndexT trackId_omvg = *s_i;
      const IndexT trackId_slamPP = slam_pp_data.getNextFreeSlamPPId();
      // Save mapping between omvg and slamPP
      slam_pp_data.setTrackId_SlamPP(trackId_omvg,trackId_slamPP);
      // Get position of landmark
      Vec3 & l_pos_w = landmarks[trackId_omvg].X;

      switch (slam_pp_data.iOutputLandmarkType)
      {
        case 0: // euclidean (world)
        {
          slam_pp_data.slamPP_DatasetFile << "VERTEX_XYZ" 
          << " " << trackId_slamPP
          << " " << l_pos_w(0)
          << " " << l_pos_w(1)
          << " " << l_pos_w(2)
          << std::endl;
        }
        break;
        case 1: // inverse depth (reference cam)
        {
          IndexT track_owner_camId_omvg = slam_pp_data.owner_track_cam_id[trackId_omvg];        
          IndexT track_owner_camId_slamPP;
          if (!slam_pp_data.getCamId_SlamPP(track_owner_camId_omvg,track_owner_camId_slamPP))
          {
            std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
          }

          const View * view_owner = sfm_data_.GetViews().at(track_owner_camId_omvg).get();
          if ( !sfm_data_.IsPoseAndIntrinsicDefined( view_owner ) )
          {
            continue;
          }
          const Pose3 pose_owner= sfm_data_.GetPoseOrDie(view_owner);
          const Mat3 rotation_owner = pose_owner.rotation();
          
          // Point to owner camera coordinate system
          Vec3 pt_owner = pose_owner(l_pos_w);
          pt_owner(0) = pt_owner(0) / pt_owner(2);
          pt_owner(1) = pt_owner(1) / pt_owner(2);
          pt_owner(2) = 1.0 / pt_owner(2);

          slam_pp_data.slamPP_DatasetFile << "VERTEX:INVD" 
          << " " << trackId_slamPP
          << " " << track_owner_camId_slamPP
          << " " << pt_owner(0)
          << " " << pt_owner(1)
          << " " << pt_owner(2)
          << std::endl;
        }
        break;
      }
    }

    // Observations - sort them by cam and track id (of slamPP)
    Hash_Map<IndexT,std::set<IndexT> >observations_slamPP;
    
    for (Hash_Map<IndexT, std::set<IndexT> >::iterator it_views_tracks_obs = increment_observation_.first.begin(); it_views_tracks_obs != increment_observation_.first.end(); ++it_views_tracks_obs)
    {
      IndexT camId_slamPP, trackId_slamPP, camId_omvg, trackId_omvg;
      // Get SlamPP index of camera
      camId_omvg = it_views_tracks_obs->first;
      if (!slam_pp_data.getCamId_SlamPP(camId_omvg,camId_slamPP))
      {
        std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
      }
      // Get iterator to the set of tracks for current view
      for (std::set<IndexT>::iterator it_view_tracks = it_views_tracks_obs->second.begin() ; it_view_tracks != it_views_tracks_obs->second.end(); ++it_view_tracks)
      {
        // Get SlamPP index of track
        trackId_omvg = *it_view_tracks;
        if (!slam_pp_data.getTrackId_SlamPP(trackId_omvg,trackId_slamPP))
        {
          std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
        }

        // observations of this cam already exist -> we just add new
        observations_slamPP[camId_slamPP].insert(trackId_slamPP);
        
      }
    }

    std::ostringstream edge_other_stream;
    slam_pp_data.initBreadthSearchFirstGraph();
    // Saved path of predecesors so we dont have to query it for every observation
    std::vector<IndexT> cam_predecesors;
    std::pair<IndexT,IndexT> last_queried_path(std::numeric_limits<IndexT>::max(),std::numeric_limits<IndexT>::max());

    // Export ordered observations to file
    for (Hash_Map<IndexT,std::set<IndexT> >::iterator it_obs = observations_slamPP.begin(); it_obs != observations_slamPP.end(); ++it_obs)
    { 
      IndexT camId_slamPP, trackId_slamPP, camId_omvg, trackId_omvg;
      camId_slamPP = it_obs->first;
      if (!slam_pp_data.getCamId_OpenMVG(camId_omvg,camId_slamPP))
      {
        std::cerr << "Something went wrong with CameraID SlamPP - OpenMVG index mapping\n";
      }
      // Get information of the camera -> to undistort the point
      const View * view = sfm_data_.GetViews().at(camId_omvg).get();
      std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic (nullptr);
      if (sfm_data_.GetIntrinsics().count(view->id_intrinsic))
      {
        optional_intrinsic = sfm_data_.GetIntrinsics().at(view->id_intrinsic);
      }

      // Loop throught measurements of the camera
      std::set<IndexT> & obs_tracks = it_obs->second;
      for (std::set<IndexT>::iterator it_tracks = obs_tracks.begin(); it_tracks != obs_tracks.end(); ++it_tracks)
      {
        trackId_slamPP = *it_tracks;
        if (!slam_pp_data.getTrackId_OpenMVG(trackId_omvg,trackId_slamPP))
        {
          std::cerr << "Something went wrong with TrackID SlamPP - OpenMVG index mapping\n";
        }

        // Get owner camid
        IndexT track_owner_camId_omvg = slam_pp_data.owner_track_cam_id[trackId_omvg];
        IndexT track_owner_camId_slamPP;
        if (!slam_pp_data.getCamId_SlamPP(track_owner_camId_omvg,track_owner_camId_slamPP))
        {
          std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
        }

        // Get feature id of the point in camId and trackId
        const IndexT featId_omvg = map_tracks_[trackId_omvg][camId_omvg];
        Vec2 pt_2d;
        pt_2d = features_provider_->feats_per_view.at(camId_omvg)[featId_omvg].coords().cast<double>();
        if (optional_intrinsic && optional_intrinsic->have_disto())
        {
          pt_2d = optional_intrinsic->get_ud_pixel(pt_2d);
        }

        switch (slam_pp_data.iOutputLandmarkType)
        {
          case 0: // Euclidean (world)
            slam_pp_data.slamPP_DatasetFile << "EDGE_PROJECT_P2MC"
            << " " << trackId_slamPP
            << " " << camId_slamPP
            << " " << pt_2d(0)
            << " " << pt_2d(1)
            << " " << "1 0 1"
            << std::endl;
          break;
          case 1:
            if (camId_slamPP == track_owner_camId_slamPP)
            {
              slam_pp_data.slamPP_DatasetFile << "EDGE_PROJ_SELF"
              << " " << trackId_slamPP
              << " " << camId_slamPP
              << " " << pt_2d(0)
              << " " << pt_2d(1)
              << " " << "1 0 1"
              << std::endl;
            }
            else
            {
              // Only query if its different from the last
              if (last_queried_path.first != camId_omvg || last_queried_path.second != track_owner_camId_omvg)
              {
                cam_predecesors.clear();
                slam_pp_data.findPathBetween(camId_omvg,track_owner_camId_omvg,cam_predecesors);
                last_queried_path.first = camId_omvg;
                last_queried_path.second = track_owner_camId_omvg;
              }
              edge_other_stream << "EDGE_PROJ_OTHER"
              << " " << trackId_slamPP
              << " " << cam_predecesors.size() + 1
              << " " << camId_slamPP;
              for (auto it_pred_camId = cam_predecesors.rbegin(); it_pred_camId != cam_predecesors.rend(); ++it_pred_camId)
              {
                IndexT pred_cam_slamPP;
                if (!slam_pp_data.getCamId_SlamPP(*it_pred_camId,pred_cam_slamPP))
                {
                  std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
                }
                edge_other_stream << " " << pred_cam_slamPP;
              }
              edge_other_stream << " " << pt_2d(0)
              << " " << pt_2d(1)
              << " " << "1 0 1"
              << std::endl;
            }
          break;
        }
      }
    }
    slam_pp_data.slamPP_DatasetFile << edge_other_stream.str();


    // Export consistency marker
    slam_pp_data.slamPP_DatasetFile << "CONSISTENCY_MARKER\n";

  }


  void IncrementalSfMReconstructionEngine::ExportTwoFoldIncrementToGraphFile_SlamPP()
  {

    // Cameras
    for (std::set<IndexT>::iterator c_i = increment_camera_.first.begin(); c_i != increment_camera_.first.end(); ++c_i)
    {
      const IndexT camId_omvg = *c_i;
      const IndexT camId_slamPP = slam_pp_data.getNextFreeSlamPPId();

      // Save mapping between omvg and slamPP
      slam_pp_data.setCamId_SlamPP(camId_omvg,camId_slamPP);

      // Get information of the camera
      const View * view = sfm_data_.GetViews().at(camId_omvg).get();
      if ( !sfm_data_.IsPoseAndIntrinsicDefined( view ) )
      {
        continue;
      }

      const IntrinsicBase * cam = sfm_data_.GetIntrinsics().at(view->id_intrinsic).get();
      const Pinhole_Intrinsic * pinhole_cam = static_cast<const Pinhole_Intrinsic *>( cam );
      const Mat3 K = pinhole_cam->K();

      switch (slam_pp_data.iOutputVertexType)
      {
        case 0: // SE(3)
        {
          // Get pose
          Pose3 pose= sfm_data_.GetPoseOrDie(view);
          const Mat3 rotation = pose.rotation().transpose();
          const Vec3 center = pose.center();
          Eigen::Quaterniond q( rotation ) ;
          // Export to graph file
          slam_pp_data.slamPP_DatasetFile << "VERTEX_CAM" 
            << " " << camId_slamPP
            << " " << center[0]
            << " " << center[1]
            << " " << center[2]
            << " " << q.x()
            << " " << q.y()
            << " " << q.z()
            << " " << q.w()
            << " " << K(0,0)
            << " " << K(1,1)
            << " " << K(0,2)
            << " " << K(1,2)
            << " " << "0.0"
            << std::endl;
        }
        break;
        case 1: // Sim(3)
        {
          Pose3 pose = sfm_data_.GetPoseOrDie(view);
          const Mat3 rotation = pose.rotation().transpose();
          const Vec3 center = pose.center();
          const double scale = 1.0;
          double angleAxis[3];
          ceres::RotationMatrixToAngleAxis((const double*)rotation.data(), angleAxis);
          // Export to graph file
          slam_pp_data.slamPP_DatasetFile << "VERTEX_CAM:SIM3" 
            << " " << camId_slamPP
            << " " << center[0]
            << " " << center[1]
            << " " << center[2]
            << " " << angleAxis[0]
            << " " << angleAxis[1]
            << " " << angleAxis[2]
            << " " << scale
            << " " << K(0,0)
            << " " << K(1,1)
            << " " << K(0,2)
            << " " << K(1,2)
            << " " << "0.0"
            << std::endl;
        }
        break;
      }
    
    }


    // Structure
    std::ostringstream new_landmarks_stream;

    Landmarks & landmarks = sfm_data_.structure;
    for (std::set<IndexT>::iterator s_i = increment_structure_.first.begin(); s_i != increment_structure_.first.end(); ++s_i)
    {
      const IndexT trackId_omvg = *s_i;
      const IndexT trackId_slamPP = slam_pp_data.getNextFreeSlamPPId();
      // Save mapping between omvg and slamPP
      slam_pp_data.setTrackId_SlamPP(trackId_omvg,trackId_slamPP);
      // Get position of landmark
      Vec3 & l_pos_w = landmarks[trackId_omvg].X;

      switch (slam_pp_data.iOutputLandmarkType)
      {
        case 0: // euclidean (world)
        {
          new_landmarks_stream << "VERTEX_XYZ" 
          << " " << trackId_slamPP
          << " " << l_pos_w(0)
          << " " << l_pos_w(1)
          << " " << l_pos_w(2)
          << std::endl;
        }
        break;
        case 1: // inverse depth (reference cam)
        {
          IndexT track_owner_camId_omvg = slam_pp_data.owner_track_cam_id[trackId_omvg];        
          IndexT track_owner_camId_slamPP;
          if (!slam_pp_data.getCamId_SlamPP(track_owner_camId_omvg,track_owner_camId_slamPP))
          {
            std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
          }

          const View * view_owner = sfm_data_.GetViews().at(track_owner_camId_omvg).get();
          if ( !sfm_data_.IsPoseAndIntrinsicDefined( view_owner ) )
          {
            continue;
          }
          const Pose3 pose_owner= sfm_data_.GetPoseOrDie(view_owner);
          const Mat3 rotation_owner = pose_owner.rotation();
          
          // Point to owner camera coordinate system
          Vec3 pt_owner = pose_owner(l_pos_w);
          pt_owner(0) = pt_owner(0) / pt_owner(2);
          pt_owner(1) = pt_owner(1) / pt_owner(2);
          pt_owner(2) = 1.0 / pt_owner(2);

          new_landmarks_stream << "VERTEX:INVD" 
          << " " << trackId_slamPP
          << " " << track_owner_camId_slamPP
          << " " << pt_owner(0)
          << " " << pt_owner(1)
          << " " << pt_owner(2)
          << std::endl;
        }
        break;
      }
    }


    slam_pp_data.initBreadthSearchFirstGraph();
    // Saved path of predecesors so we dont have to query it for every observation
    std::vector<IndexT> cam_predecesors;
    std::pair<IndexT,IndexT> last_queried_path(std::numeric_limits<IndexT>::max(),std::numeric_limits<IndexT>::max());

    // Observations (of old points) - sort them by cam and track id (of slamPP)
    Hash_Map<IndexT,std::set<IndexT> >observations_slamPP;
    
    for (Hash_Map<IndexT, std::set<IndexT> >::iterator it_views_tracks_obs = increment_observation_.first.begin(); it_views_tracks_obs != increment_observation_.first.end(); ++it_views_tracks_obs)
    {
      IndexT camId_slamPP, trackId_slamPP, camId_omvg, trackId_omvg;
      // Get SlamPP index of camera
      camId_omvg = it_views_tracks_obs->first;
      if (!slam_pp_data.getCamId_SlamPP(camId_omvg,camId_slamPP))
      {
        std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
      }
      // Get iterator to the set of tracks for current view
      for (std::set<IndexT>::iterator it_view_tracks = it_views_tracks_obs->second.begin() ; it_view_tracks != it_views_tracks_obs->second.end(); ++it_view_tracks)
      {
        // Get SlamPP index of track
        trackId_omvg = *it_view_tracks;
        if (!slam_pp_data.getTrackId_SlamPP(trackId_omvg,trackId_slamPP))
        {
          std::cerr << "Something went wrong with Track ID OpenMVG B- SlamPP index mapping\n";
        }

        // Skip if its a new landmarks - we only add observations of already reconstructed points
        if (increment_structure_.first.find(trackId_omvg) != increment_structure_.first.end())
        {
          continue;
        }

        // observations of this cam already exist -> we just add new
        observations_slamPP[camId_slamPP].insert(trackId_slamPP);
        
      }
    }

    // Observations of the points already reconstructed (without the newly added camera)

    // Export ordered observations to file
    for (Hash_Map<IndexT,std::set<IndexT> >::iterator it_obs = observations_slamPP.begin(); it_obs != observations_slamPP.end(); ++it_obs)
    { 
      IndexT camId_slamPP, trackId_slamPP, camId_omvg, trackId_omvg;
      camId_slamPP = it_obs->first;
      if (!slam_pp_data.getCamId_OpenMVG(camId_omvg,camId_slamPP))
      {
        std::cerr << "Something went wrong with CameraID SlamPP - OpenMVG index mapping\n";
      }
      // Get information of the camera -> to undistort the point
      const View * view = sfm_data_.GetViews().at(camId_omvg).get();
      std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic (nullptr);
      if (sfm_data_.GetIntrinsics().count(view->id_intrinsic))
      {
        optional_intrinsic = sfm_data_.GetIntrinsics().at(view->id_intrinsic);
      }

      // Loop throught measurements of the camera
      std::set<IndexT> & obs_tracks = it_obs->second;
      for (std::set<IndexT>::iterator it_tracks = obs_tracks.begin(); it_tracks != obs_tracks.end(); ++it_tracks)
      {
        trackId_slamPP = *it_tracks;
        if (!slam_pp_data.getTrackId_OpenMVG(trackId_omvg,trackId_slamPP))
        {
          std::cerr << "Something went wrong with TrackID SlamPP - OpenMVG index mapping\n";
        }

        // Get owner camid
        IndexT track_owner_camId_omvg = slam_pp_data.owner_track_cam_id[trackId_omvg];
        IndexT track_owner_camId_slamPP;
        if (!slam_pp_data.getCamId_SlamPP(track_owner_camId_omvg,track_owner_camId_slamPP))
        {
          std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
        }

        // Get feature id of the point in camId and trackId
        const IndexT featId_omvg = map_tracks_[trackId_omvg][camId_omvg];
        Vec2 pt_2d;
        pt_2d = features_provider_->feats_per_view.at(camId_omvg)[featId_omvg].coords().cast<double>();
        if (optional_intrinsic && optional_intrinsic->have_disto())
        {
          pt_2d = optional_intrinsic->get_ud_pixel(pt_2d);
        }

        switch (slam_pp_data.iOutputLandmarkType)
        {
          case 0: // Euclidean (world)
            slam_pp_data.slamPP_DatasetFile << "EDGE_PROJECT_P2MC"
            << " " << trackId_slamPP
            << " " << camId_slamPP
            << " " << pt_2d(0)
            << " " << pt_2d(1)
            << " " << "1 0 1"
            << std::endl;
          break;
          case 1:
            if (camId_slamPP == track_owner_camId_slamPP)
            {
              // Actually shouldnt happen because the point should already be in
              slam_pp_data.slamPP_DatasetFile << "EDGE_PROJ_SELF"
              << " " << trackId_slamPP
              << " " << camId_slamPP
              << " " << pt_2d(0)
              << " " << pt_2d(1)
              << " " << "1 0 1"
              << std::endl;
            }
            else
            {
              // Only query if its different from the last
              if (last_queried_path.first != camId_omvg || last_queried_path.second != track_owner_camId_omvg)
              {
                cam_predecesors.clear();
                slam_pp_data.findPathBetween(camId_omvg,track_owner_camId_omvg,cam_predecesors);
                last_queried_path.first = camId_omvg;
                last_queried_path.second = track_owner_camId_omvg;
              }
              slam_pp_data.slamPP_DatasetFile << "EDGE_PROJ_OTHER"
              << " " << trackId_slamPP
              << " " << cam_predecesors.size() + 1
              << " " << camId_slamPP;
              for (auto it_pred_camId = cam_predecesors.rbegin(); it_pred_camId != cam_predecesors.rend(); ++it_pred_camId)
              {
                IndexT pred_cam_slamPP;
                if (!slam_pp_data.getCamId_SlamPP(*it_pred_camId,pred_cam_slamPP))
                {
                  std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
                }
                slam_pp_data.slamPP_DatasetFile << " " << pred_cam_slamPP;
              }
              slam_pp_data.slamPP_DatasetFile << " " << pt_2d(0)
              << " " << pt_2d(1)
              << " " << "1 0 1"
              << std::endl;
            }
          break;
        }
      }
    }


    // Export consistency marker if there were any observations of old points added
    if (observations_slamPP.size()!=0)
    {
      slam_pp_data.slamPP_DatasetFile << "CONSISTENCY_MARKER\n";
    }
    
    slam_pp_data.slamPP_DatasetFile << new_landmarks_stream.str();

    // Add the observations of new points in reference image and other observations of new points
    observations_slamPP.clear();
    
    for (Hash_Map<IndexT, std::set<IndexT> >::iterator it_views_tracks_obs = increment_observation_.first.begin(); it_views_tracks_obs != increment_observation_.first.end(); ++it_views_tracks_obs)
    {
      IndexT camId_slamPP, trackId_slamPP, camId_omvg, trackId_omvg;
      // Get SlamPP index of camera
      camId_omvg = it_views_tracks_obs->first;

      if (!slam_pp_data.getCamId_SlamPP(camId_omvg,camId_slamPP))
      {
        std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
      }
      // Get iterator to the set of tracks for current view
      for (std::set<IndexT>::iterator it_view_tracks = it_views_tracks_obs->second.begin() ; it_view_tracks != it_views_tracks_obs->second.end(); ++it_view_tracks)
      {
        // Get SlamPP index of track
        trackId_omvg = *it_view_tracks;
        if (!slam_pp_data.getTrackId_SlamPP(trackId_omvg,trackId_slamPP))
        {
          std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
        }

        // Skip if its an old landmark (already added)
        if (increment_structure_.first.find(trackId_omvg) == increment_structure_.first.end())
        {
          continue;
        }

        // observations of this cam already exist -> we just add new
        observations_slamPP[camId_slamPP].insert(trackId_slamPP);
      }
    }

    std::ostringstream edge_other_stream;

    // Export ordered observations to file
    for (Hash_Map<IndexT,std::set<IndexT> >::iterator it_obs = observations_slamPP.begin(); it_obs != observations_slamPP.end(); ++it_obs)
    { 
      IndexT camId_slamPP, trackId_slamPP, camId_omvg, trackId_omvg;
      camId_slamPP = it_obs->first;
      if (!slam_pp_data.getCamId_OpenMVG(camId_omvg,camId_slamPP))
      {
        std::cerr << "Something went wrong with CameraID SlamPP - OpenMVG index mapping\n";
      }
      // Get information of the camera -> to undistort the point
      const View * view = sfm_data_.GetViews().at(camId_omvg).get();
      std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic (nullptr);
      if (sfm_data_.GetIntrinsics().count(view->id_intrinsic))
      {
        optional_intrinsic = sfm_data_.GetIntrinsics().at(view->id_intrinsic);
      }

      // Loop throught measurements of the camera
      std::set<IndexT> & obs_tracks = it_obs->second;
      for (std::set<IndexT>::iterator it_tracks = obs_tracks.begin(); it_tracks != obs_tracks.end(); ++it_tracks)
      {
        trackId_slamPP = *it_tracks;
        if (!slam_pp_data.getTrackId_OpenMVG(trackId_omvg,trackId_slamPP))
        {
          std::cerr << "Something went wrong with TrackID SlamPP - OpenMVG index mapping\n";
        }

        // Get owner camid
        IndexT track_owner_camId_omvg = slam_pp_data.owner_track_cam_id[trackId_omvg];
        IndexT track_owner_camId_slamPP;
        if (!slam_pp_data.getCamId_SlamPP(track_owner_camId_omvg,track_owner_camId_slamPP))
        {
          std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
        }

        // Get feature id of the point in camId and trackId
        const IndexT featId_omvg = map_tracks_[trackId_omvg][camId_omvg];
        Vec2 pt_2d;
        pt_2d = features_provider_->feats_per_view.at(camId_omvg)[featId_omvg].coords().cast<double>();
        if (optional_intrinsic && optional_intrinsic->have_disto())
        {
          pt_2d = optional_intrinsic->get_ud_pixel(pt_2d);
        }

        switch (slam_pp_data.iOutputLandmarkType)
        {
          case 0: // Euclidean (world)
            slam_pp_data.slamPP_DatasetFile << "EDGE_PROJECT_P2MC"
            << " " << trackId_slamPP
            << " " << camId_slamPP
            << " " << pt_2d(0)
            << " " << pt_2d(1)
            << " " << "1 0 1"
            << std::endl;
          break;
          case 1:
            if (camId_slamPP == track_owner_camId_slamPP)
            {
              slam_pp_data.slamPP_DatasetFile << "EDGE_PROJ_SELF"
              << " " << trackId_slamPP
              << " " << camId_slamPP
              << " " << pt_2d(0)
              << " " << pt_2d(1)
              << " " << "1 0 1"
              << std::endl;
            }
            else
            {
              // Only query if its different from the last
              if (last_queried_path.first != camId_omvg || last_queried_path.second != track_owner_camId_omvg)
              {
                cam_predecesors.clear();
                slam_pp_data.findPathBetween(camId_omvg,track_owner_camId_omvg,cam_predecesors);
                last_queried_path.first = camId_omvg;
                last_queried_path.second = track_owner_camId_omvg;
              }
              edge_other_stream << "EDGE_PROJ_OTHER"
              << " " << trackId_slamPP
              << " " << cam_predecesors.size() + 1
              << " " << camId_slamPP;
              for (auto it_pred_camId = cam_predecesors.rbegin(); it_pred_camId != cam_predecesors.rend(); ++it_pred_camId)
              {
                IndexT pred_cam_slamPP;
                if (!slam_pp_data.getCamId_SlamPP(*it_pred_camId,pred_cam_slamPP))
                {
                  std::cerr << "Something went wrong with Camera ID OpenMVG - SlamPP index mapping\n";
                }
                edge_other_stream << " " << pred_cam_slamPP;
              }
              edge_other_stream << " " << pt_2d(0)
              << " " << pt_2d(1)
              << " " << "1 0 1"
              << std::endl;
            }
          break;
        }
      }
    }
    slam_pp_data.slamPP_DatasetFile << edge_other_stream.str();


    // Export consistency marker
    slam_pp_data.slamPP_DatasetFile << "CONSISTENCY_MARKER\n";

  }

*/

  void IncrementalSfMReconstructionEngine::ExportTwoFoldTotalProcessByIncrementToGraphFile_SlamPP()
  {

    // Create output files
    std::string graphfile = stlplus::create_filespec(slam_pp_data.graphOutputDir, "GraphFile", ".txt");
    std::string camfile = stlplus::create_filespec(slam_pp_data.graphOutputDir, "CAM_OpenMVG_SlamPP", ".txt");
    std::string pointsfile = stlplus::create_filespec(slam_pp_data.graphOutputDir, "3D_OpenMVG_SlamPP", ".txt");

    // Open graph file
    std::ofstream slamPP_GraphFile;
    slamPP_GraphFile.open( graphfile.c_str(), std::ios::out );

    SlamPP_Data slam_pp_total_data;
    slam_pp_total_data.initCamParentsGraph();

    for (IndexT i_inc_iter = 0; i_inc_iter < history_i_camera_.size(); i_inc_iter++)
    {
      std::set<IndexT> & current_inc_camera = history_i_camera_[i_inc_iter].first;
      std::set<IndexT> & current_inc_structure = history_i_structure_[i_inc_iter].first;
      Hash_Map<IndexT, std::set<IndexT> > & current_inc_observation = history_i_observations_[i_inc_iter].first;


      // Cameras
      for (std::set<IndexT>::iterator c_i = current_inc_camera.begin(); c_i != current_inc_camera.end(); ++c_i)
      {
        const IndexT camId_omvg = *c_i;
        const IndexT camId_slamPP = slam_pp_total_data.getNextFreeSlamPPId();

        // Add camera to parents graph
        if (i_inc_iter !=0)
          slam_pp_total_data.addCamWithParentToGraph(camId_omvg, slam_pp_data.parent_cam_id[camId_omvg]);

        // Save mapping between omvg and slamPP
        slam_pp_total_data.setCamId_SlamPP(camId_omvg,camId_slamPP);

        // Get information of the camera
        const View * view = sfm_data_.GetViews().at(camId_omvg).get();
        if ( !sfm_data_.IsPoseAndIntrinsicDefined( view ) )
        {
          continue;
        }

        const IntrinsicBase * cam = sfm_data_.GetIntrinsics().at(view->id_intrinsic).get();
        const Pinhole_Intrinsic * pinhole_cam = static_cast<const Pinhole_Intrinsic *>( cam );
        const Mat3 K = pinhole_cam->K();

        switch (slam_pp_data.iOutputVertexType)
        {
          case 0: // SE(3)
          {
            // Get pose
            Pose3 pose= sfm_data_.GetPoseOrDie(view);
            const Mat3 rotation = pose.rotation().transpose();
            const Vec3 center = pose.center();
            Eigen::Quaterniond q( rotation ) ;
            // Export to graph file
            slamPP_GraphFile << "VERTEX_CAM" 
              << " " << camId_slamPP
              << " " << center[0]
              << " " << center[1]
              << " " << center[2]
              << " " << q.x()
              << " " << q.y()
              << " " << q.z()
              << " " << q.w()
              << " " << K(0,0)
              << " " << K(1,1)
              << " " << K(0,2)
              << " " << K(1,2)
              << " " << "0.0"
              << std::endl;
          }
          break;
          case 1: // Sim(3)
          {
            Pose3 pose = sfm_data_.GetPoseOrDie(view);
            const Mat3 rotation = pose.rotation().transpose();
            const Vec3 center = pose.center();
            const double scale = 1.0;
            double angleAxis[3];
            ceres::RotationMatrixToAngleAxis((const double*)rotation.data(), angleAxis);
            // Export to graph file
            slamPP_GraphFile << "VERTEX_CAM:SIM3" 
              << " " << camId_slamPP
              << " " << center[0]
              << " " << center[1]
              << " " << center[2]
              << " " << angleAxis[0]
              << " " << angleAxis[1]
              << " " << angleAxis[2]
              << " " << scale
              << " " << K(0,0)
              << " " << K(1,1)
              << " " << K(0,2)
              << " " << K(1,2)
              << " " << "0.0"
              << std::endl;
          }
          break;
        }
      }

      if(i_inc_iter == 0)
      {
        for (std::set<IndexT>::reverse_iterator c_i = current_inc_camera.rbegin(); c_i != current_inc_camera.rend(); ++c_i)
        {
            // Get SlamPP index of camera
            IndexT camId_omvg = *c_i;
            slam_pp_total_data.addCamWithParentToGraph(camId_omvg, slam_pp_data.parent_cam_id[camId_omvg]);
        }
      }
      // Structure
      std::ostringstream new_landmarks_stream;

      Landmarks & landmarks = sfm_data_.structure;
      for (std::set<IndexT>::iterator s_i = current_inc_structure.begin(); s_i != current_inc_structure.end(); ++s_i)
      {
        const IndexT trackId_omvg = *s_i;

        // If point was ever removed and the removal was in later iteration that currently inspected --> we skip the structure as it will be removed
        if (structure_last_removed.count(trackId_omvg) != 0 && !(i_inc_iter > structure_last_removed[trackId_omvg]))
          continue;

        const IndexT trackId_slamPP = slam_pp_total_data.getNextFreeSlamPPId();
        // Save mapping between omvg and slamPP
        slam_pp_total_data.setTrackId_SlamPP(trackId_omvg,trackId_slamPP);
        // Get position of landmark
        Vec3 & l_pos_w = landmarks[trackId_omvg].X;

        switch (slam_pp_data.iOutputLandmarkType)
        {
          case 0: // euclidean (world)
          {
            new_landmarks_stream << "VERTEX_XYZ" 
            << " " << trackId_slamPP
            << " " << l_pos_w(0)
            << " " << l_pos_w(1)
            << " " << l_pos_w(2)
            << std::endl;
          }
          break;
          case 1: // inverse depth (reference cam)
          {
            IndexT track_owner_camId_omvg = slam_pp_data.owner_track_cam_id[trackId_omvg];        
            IndexT track_owner_camId_slamPP;
            if (!slam_pp_total_data.getCamId_SlamPP(track_owner_camId_omvg,track_owner_camId_slamPP))
            {
              std::cerr << "Something went wrong with A Camera ID OpenMVG - SlamPP index mapping "<<track_owner_camId_omvg<<"\n";
            }

            const View * view_owner = sfm_data_.GetViews().at(track_owner_camId_omvg).get();
            if ( !sfm_data_.IsPoseAndIntrinsicDefined( view_owner ) )
            {
              continue;
            }
            const Pose3 pose_owner= sfm_data_.GetPoseOrDie(view_owner);
            const Mat3 rotation_owner = pose_owner.rotation();
            
            // Point to owner camera coordinate system
            Vec3 pt_owner = pose_owner(l_pos_w);
            pt_owner(0) = pt_owner(0) / pt_owner(2);
            pt_owner(1) = pt_owner(1) / pt_owner(2);
            pt_owner(2) = 1.0 / pt_owner(2);

            new_landmarks_stream << "VERTEX:INVD" 
            << " " << trackId_slamPP
            << " " << track_owner_camId_slamPP
            << " " << pt_owner(0)
            << " " << pt_owner(1)
            << " " << pt_owner(2)
            << std::endl;
          }
          break;
        }
      }

      slam_pp_total_data.initBreadthSearchFirstGraph();
      // Saved path of predecesors so we dont have to query it for every observation
      std::vector<IndexT> cam_predecesors;
      std::pair<IndexT,IndexT> last_queried_path(std::numeric_limits<IndexT>::max(),std::numeric_limits<IndexT>::max());

      // Observations (of old points) - sort them by cam and track id (of slamPP)
      Hash_Map<IndexT,std::set<IndexT> > oldPoints_observations_slamPP;
      // Observations (of new points) - sort them by cam and track id (of slamPP)
      Hash_Map<IndexT,std::set<IndexT> > newPoints_observations_slamPP;

      for (Hash_Map<IndexT, std::set<IndexT> >::iterator it_views_tracks_obs = current_inc_observation.begin(); it_views_tracks_obs != current_inc_observation.end(); ++it_views_tracks_obs)
      {
        IndexT camId_slamPP, trackId_slamPP, camId_omvg, trackId_omvg;
        // Get SlamPP index of camera
        camId_omvg = it_views_tracks_obs->first;
        if (!slam_pp_total_data.getCamId_SlamPP(camId_omvg,camId_slamPP))
        {
          std::cerr << "Something went wrong with B Camera ID OpenMVG - SlamPP index mapping"<<camId_omvg<<"\n";
        }
        // Get iterator to the set of tracks for current view
        for (std::set<IndexT>::iterator it_view_tracks = it_views_tracks_obs->second.begin() ; it_view_tracks != it_views_tracks_obs->second.end(); ++it_view_tracks)
        {
          // Get SlamPP index of track
          trackId_omvg = *it_view_tracks;

          // If structue was ever removed and the removal was in later iteration that currently inspected --> we skip the structure as it will be removed
          if (structure_last_removed.count(trackId_omvg) != 0 && !(i_inc_iter > structure_last_removed[trackId_omvg]))
            continue;
      
          // Check if this measurement was deleted after this iteration
          if (observation_last_removed.count(camId_omvg) != 0 && observation_last_removed[camId_omvg].count(trackId_omvg) != 0 && !(i_inc_iter > observation_last_removed[camId_omvg][trackId_omvg]))
            continue;    

          if (!slam_pp_total_data.getTrackId_SlamPP(trackId_omvg,trackId_slamPP))
          {
            std::cerr << "Something went wrong with C Track ID OpenMVG A - SlamPP index mapping"<<trackId_omvg<<"\n";
          }
          else
          {
            // Skip if its a new landmarks - we only add observations of already reconstructed points
            if (current_inc_structure.find(trackId_omvg) == current_inc_structure.end())
            {
              // observations of this cam already exist -> we just add new
              oldPoints_observations_slamPP[camId_slamPP].insert(trackId_slamPP);
            }
            else
            {
              // observations of this cam already exist -> we just add new
              newPoints_observations_slamPP[camId_slamPP].insert(trackId_slamPP);
            }
          }
        }
      }

      // Export ordered observations of old points to file
      for (Hash_Map<IndexT,std::set<IndexT> >::iterator it_obs = oldPoints_observations_slamPP.begin(); it_obs != oldPoints_observations_slamPP.end(); ++it_obs)
      { 
        IndexT camId_slamPP, trackId_slamPP, camId_omvg, trackId_omvg;
        camId_slamPP = it_obs->first;
        if (!slam_pp_total_data.getCamId_OpenMVG(camId_omvg,camId_slamPP))
        {
          std::cerr << "Something went wrong with E CameraID SlamPP - OpenMVG index mapping"<<camId_slamPP<<"\n";
        }
        // Get information of the camera -> to undistort the point
        const View * view = sfm_data_.GetViews().at(camId_omvg).get();
        std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic (nullptr);
        if (sfm_data_.GetIntrinsics().count(view->id_intrinsic))
        {
          optional_intrinsic = sfm_data_.GetIntrinsics().at(view->id_intrinsic);
        }

        // Loop throught measurements of the camera
        std::set<IndexT> & obs_tracks = it_obs->second;
        for (std::set<IndexT>::iterator it_tracks = obs_tracks.begin(); it_tracks != obs_tracks.end(); ++it_tracks)
        {
          trackId_slamPP = *it_tracks;
          if (!slam_pp_total_data.getTrackId_OpenMVG(trackId_omvg,trackId_slamPP))
          {
            std::cerr << "Something went wrong with F TrackID SlamPP - OpenMVG index mapping"<<trackId_slamPP<<"\n";
          }

          // Get owner camid
          IndexT track_owner_camId_omvg = slam_pp_data.owner_track_cam_id[trackId_omvg];
          IndexT track_owner_camId_slamPP;
          if (!slam_pp_total_data.getCamId_SlamPP(track_owner_camId_omvg,track_owner_camId_slamPP))
          {
            std::cerr << "Something went wrong with G Camera ID OpenMVG - SlamPP index mapping"<<track_owner_camId_omvg<<"\n";
          }

          // Get feature id of the point in camId and trackId
          const IndexT featId_omvg = map_tracks_[trackId_omvg][camId_omvg];
          Vec2 pt_2d;
          pt_2d = features_provider_->feats_per_view.at(camId_omvg)[featId_omvg].coords().cast<double>();
          if (optional_intrinsic && optional_intrinsic->have_disto())
          {
            pt_2d = optional_intrinsic->get_ud_pixel(pt_2d);
          }

          switch (slam_pp_data.iOutputLandmarkType)
          {
            case 0: // Euclidean (world)
              slamPP_GraphFile << "EDGE_PROJECT_P2MC"
              << " " << trackId_slamPP
              << " " << camId_slamPP
              << " " << pt_2d(0)
              << " " << pt_2d(1)
              << " " << "1 0 1"
              << std::endl;
            break;
            case 1:
              if (camId_slamPP == track_owner_camId_slamPP)
              {
                // Actually shouldnt happen because the point should already be in
                slamPP_GraphFile << "EDGE_PROJ_SELF"
                << " " << trackId_slamPP
                << " " << camId_slamPP
                << " " << pt_2d(0)
                << " " << pt_2d(1)
                << " " << "1 0 1"
                << std::endl;
              }
              else
              {
                // Only query if its different from the last
                if (last_queried_path.first != camId_omvg || last_queried_path.second != track_owner_camId_omvg)
                {
                  cam_predecesors.clear();
                  slam_pp_total_data.findPathBetween(camId_omvg,track_owner_camId_omvg,cam_predecesors);
                  last_queried_path.first = camId_omvg;
                  last_queried_path.second = track_owner_camId_omvg;
                }
                slamPP_GraphFile << "EDGE_PROJ_OTHER"
                << " " << trackId_slamPP
                << " " << cam_predecesors.size() + 1
                << " " << camId_slamPP;
                for (auto it_pred_camId = cam_predecesors.rbegin(); it_pred_camId != cam_predecesors.rend(); ++it_pred_camId)
                {
                  IndexT pred_cam_slamPP;
                  if (!slam_pp_total_data.getCamId_SlamPP(*it_pred_camId,pred_cam_slamPP))
                  {
                    std::cerr << "Something went wrong with H Camera ID OpenMVG - SlamPP index mapping"<<*it_pred_camId<<"\n";
                  }
                  slamPP_GraphFile << " " << pred_cam_slamPP;
                }
                slamPP_GraphFile << " " << pt_2d(0)
                << " " << pt_2d(1)
                << " " << "1 0 1"
                << std::endl;
              }
            break;
          }
        }
      }

      // Export consistency marker if there were any observations of old points added
      if (oldPoints_observations_slamPP.size()!=0)
      {
        slamPP_GraphFile << "CONSISTENCY_MARKER\n";
      }
      
      // Dump new stucture points to the graph file
      slamPP_GraphFile << new_landmarks_stream.str();

      // First we add the measurements on the owner camera and then the other measurements of the point ("EDGE_OTHER")
      std::ostringstream edge_other_stream;

      // Export ordered observations of new points to file
      for (Hash_Map<IndexT,std::set<IndexT> >::iterator it_obs = newPoints_observations_slamPP.begin(); it_obs != newPoints_observations_slamPP.end(); ++it_obs)
      { 
        IndexT camId_slamPP, trackId_slamPP, camId_omvg, trackId_omvg;
        camId_slamPP = it_obs->first;
        if (!slam_pp_total_data.getCamId_OpenMVG(camId_omvg,camId_slamPP))
        {
          std::cerr << "Something went wrong with K CameraID SlamPP - OpenMVG index mapping"<<camId_slamPP<<"\n";
        }
        // Get information of the camera -> to undistort the point
        const View * view = sfm_data_.GetViews().at(camId_omvg).get();
        std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic (nullptr);
        if (sfm_data_.GetIntrinsics().count(view->id_intrinsic))
        {
          optional_intrinsic = sfm_data_.GetIntrinsics().at(view->id_intrinsic);
        }

        // Loop throught measurements of the camera
        std::set<IndexT> & obs_tracks = it_obs->second;
        for (std::set<IndexT>::iterator it_tracks = obs_tracks.begin(); it_tracks != obs_tracks.end(); ++it_tracks)
        {
          trackId_slamPP = *it_tracks;
          if (!slam_pp_total_data.getTrackId_OpenMVG(trackId_omvg,trackId_slamPP))
          {
            std::cerr << "Something went wrong with J TrackID SlamPP - OpenMVG index mapping"<<trackId_slamPP<<"\n";
          }

          // Get owner camid
          IndexT track_owner_camId_omvg = slam_pp_data.owner_track_cam_id[trackId_omvg];
          IndexT track_owner_camId_slamPP;
          if (!slam_pp_total_data.getCamId_SlamPP(track_owner_camId_omvg,track_owner_camId_slamPP))
          {
            std::cerr << "Something went wrong with K Camera ID OpenMVG - SlamPP index mapping"<<track_owner_camId_omvg<<"\n";
          }

          // Get feature id of the point in camId and trackId
          const IndexT featId_omvg = map_tracks_[trackId_omvg][camId_omvg];
          Vec2 pt_2d;
          pt_2d = features_provider_->feats_per_view.at(camId_omvg)[featId_omvg].coords().cast<double>();
          if (optional_intrinsic && optional_intrinsic->have_disto())
          {
            pt_2d = optional_intrinsic->get_ud_pixel(pt_2d);
          }

          switch (slam_pp_data.iOutputLandmarkType)
          {
            case 0: // Euclidean (world)
              slamPP_GraphFile << "EDGE_PROJECT_P2MC"
              << " " << trackId_slamPP
              << " " << camId_slamPP
              << " " << pt_2d(0)
              << " " << pt_2d(1)
              << " " << "1 0 1"
              << std::endl;
            break;
            case 1:
              if (camId_slamPP == track_owner_camId_slamPP)
              {
                slamPP_GraphFile << "EDGE_PROJ_SELF"
                << " " << trackId_slamPP
                << " " << camId_slamPP
                << " " << pt_2d(0)
                << " " << pt_2d(1)
                << " " << "1 0 1"
                << std::endl;
              }
              else
              {
                // Only query if its different from the last
                if (last_queried_path.first != camId_omvg || last_queried_path.second != track_owner_camId_omvg)
                {
                  cam_predecesors.clear();
                  slam_pp_total_data.findPathBetween(camId_omvg,track_owner_camId_omvg,cam_predecesors);
                  last_queried_path.first = camId_omvg;
                  last_queried_path.second = track_owner_camId_omvg;
                }
                edge_other_stream << "EDGE_PROJ_OTHER"
                << " " << trackId_slamPP
                << " " << cam_predecesors.size() + 1
                << " " << camId_slamPP;
                for (auto it_pred_camId = cam_predecesors.rbegin(); it_pred_camId != cam_predecesors.rend(); ++it_pred_camId)
                {
                  IndexT pred_cam_slamPP;
                  if (!slam_pp_total_data.getCamId_SlamPP(*it_pred_camId,pred_cam_slamPP))
                  {
                    std::cerr << "Something went wrong with L Camera ID OpenMVG - SlamPP index mapping"<<*it_pred_camId<<"\n";
                  }
                  edge_other_stream << " " << pred_cam_slamPP;
                }
                edge_other_stream << " " << pt_2d(0)
                << " " << pt_2d(1)
                << " " << "1 0 1"
                << std::endl;
              }
            break;
          }
        }
      }
      // Add the data of measurements of non-owner cameras to the graphfile
      slamPP_GraphFile << edge_other_stream.str();

      // Export consistency marker
      slamPP_GraphFile << "CONSISTENCY_MARKER\n";
    }

    // Close the graphfile
    slamPP_GraphFile.flush();
    slamPP_GraphFile.close();


    // Export camera information to a file
    std::ofstream slamPP_CamFile;
    slamPP_CamFile.open( camfile.c_str(), std::ios::out );
    
    // Header
    slamPP_CamFile << "OpenMVG_ID,SlamPP_ID,image_filename\n";

    // Loop through cameras
    for (auto cam_data : slam_pp_total_data.camera_ids_slamPP_omvg)
    {
      const View * view = sfm_data_.GetViews().at(cam_data.second).get();
      slamPP_CamFile << cam_data.second <<","<<cam_data.first<<","<< view->s_Img_path.substr(1,view->s_Img_path.size()-1)<<"\n";
    }

    // Close the camfile
    slamPP_CamFile.flush();
    slamPP_CamFile.close();


    // Export structure information to a file
    std::ofstream slamPP_3DFile;
    slamPP_3DFile.open( pointsfile.c_str(), std::ios::out );
    
    // Header
    slamPP_3DFile << "OpenMVG_ID,SlamPP_ID\n";

    // Loop through tracks
    for (auto struct_data : slam_pp_total_data.track_ids_slamPP_omvg)
    {
      slamPP_3DFile << struct_data.second <<","<<struct_data.first<<"\n";
    }

    // Close the 3dfile
    slamPP_3DFile.flush();
    slamPP_3DFile.close();


  }



} // namespace sfm
} // namespace openMVG

