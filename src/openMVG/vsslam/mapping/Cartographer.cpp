
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include <openMVG/vsslam/mapping/Cartographer.hpp>

using namespace openMVG;

namespace openMVG  {
namespace VSSLAM  {

// Constructor
Cartographer::Cartographer(){}


void Cartographer::clearInitializationData()
{
  init_map_frames.clear();
  map_initialized = false;
}

bool Cartographer::initializationAddStep
(
  std::shared_ptr<Frame> & frame,
  std::vector<std::unique_ptr<MapLandmark> > * vec_new_pts_3D_obs
)
{
  step_id++;

  const size_t & frame_id = frame->getFrameId();

  std::cout<<"Cartographer: Add initialization step "<<frame_id<<"\n";

  if (!map_initialized)
  {
    // Add frame to the list of initial map frames
    init_map_frames.emplace_back(frame->share_ptr());

    // Add observations to local points
    // (we increase the min connectivity requirement to prevent adding to global system)
    addObservationsToLandmarks(frame.get(), min_obs_per_landmark+1);

    // Add new points to local map
    // (we increase the min connectivity requirement to prevent adding to global system)
    if (vec_new_pts_3D_obs)
      addLandmarksToStructure(frame.get(),*vec_new_pts_3D_obs, min_obs_per_landmark+1);

    // Once we have enough frames we try to initialize map
    if (init_map_frames.size() == min_obs_per_landmark)
    {
      std::cout<<"Cartographer: Initialize map\n";

      // Find min landmark degree necessary for map initialization with init_min_points
      size_t min_degree_connectivity = findMinLandmarkDegreeForGlobalMapInitialization();

      if (min_degree_connectivity==0)
      {
        std::cout<<"Cartographer: Map Initialization FAILED\n";
        clearInitializationData();
        return false;
      }

      // --------------------
      // -- Add frames to global map
      // --------------------
      for (auto & imf : init_map_frames)
      {
        addFrameToGlobalMap(imf);
      }

      // --------------------
      // -- Add landmarks which satisfy min_degree_connectivity add from local to global
      // --------------------
      for (Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> >::iterator iter_ml = tmp_structure.begin(); iter_ml!=tmp_structure.end();)
      {
        if ((iter_ml->first)->isValidByConnectivityDegree(min_degree_connectivity))
        {
          addLandmarkFromLocalToGlobalMap(iter_ml->second);
          iter_ml = tmp_structure.erase(iter_ml);
        }
        else
        {
          ++iter_ml;
        }
      }

      // Clear all unnecessary data
      clearInitializationData();
      // Set map as initialized
      map_initialized = true;
      // Show some stats
      std::cout<<"MAP: initialized: F: "<<keyframes.size()<<" P: "<<structure.size()<<" T: "<<tmp_structure.size()<<"\n";
    }
  }
  return true;
}

size_t Cartographer::findMinLandmarkDegreeForGlobalMapInitialization()
{
  // Map with number of points for each landmark degree
  std::map<size_t,size_t> map_landmarkDegree_nPts;

  // Loop through points already in local/global map
  for (auto & ml : tmp_structure)
  {
    map_landmarkDegree_nPts[ml.second->obs_.size()]++;
  }

  // Go through find the first degree of landmark connectivity that has sufficient number of points
  size_t n_totalPts = 0;
  size_t min_necessary_connectivity = 0;

  for (auto iter = map_landmarkDegree_nPts.rbegin(); iter != map_landmarkDegree_nPts.rend(); ++iter)
  {
    n_totalPts += iter->second;
    if (n_totalPts >= init_min_map_pts)
    {
      // The highest degree of connectivity we will return is the min required
      min_necessary_connectivity = std::min<size_t>(min_obs_per_landmark,(iter->first));
      break;
    }
  }

  return min_necessary_connectivity;
}

size_t Cartographer::findMinLandmarkDegreeForDefinedInGlobalMap
(
  Frame * frame,
  std::vector<std::unique_ptr<MapLandmark> > * vec_new_pts_3D_obs
)
{
  // Map with number of points for each landmark degree
  std::map<size_t,size_t> map_landmarkDegree_nPts;

  // Loop through points already in local/global map
  for (MapLandmark * & ml : frame->map_points_)
  {
    if(!ml)
      continue;

    map_landmarkDegree_nPts[ml->obs_.size()]++;
  }

  // Loop through new triangulated points
  //  - we have to check if point has measurements from this frame
  if (vec_new_pts_3D_obs)
  {
    const size_t frame_id = frame->getFrameId();
    for (std::unique_ptr<MapLandmark> & ml : *vec_new_pts_3D_obs)
    {
      if (ml->hasFrameObservation(frame_id))
        map_landmarkDegree_nPts[ml->obs_.size()]++;
    }
  }

  // Go through find the first degree of landmark connectivity that has sufficient number of points
  size_t n_totalPts = 0;
  size_t min_necessary_connectivity = 0;

  for (auto iter = map_landmarkDegree_nPts.rbegin(); iter != map_landmarkDegree_nPts.rend(); ++iter)
  {
    n_totalPts += iter->second;
    if (n_totalPts >= min_landmark_per_frame)
    {
      // The highest degree of connectivity we will return is the min required
      // (connectivity of landmarks does not include yet this frame so we must increase it for one)
      min_necessary_connectivity = std::min<size_t>(min_obs_per_landmark,(iter->first+1));
      break;
    }
  }

  return min_necessary_connectivity;
}

void Cartographer::addFrameToGlobalMap(const std::shared_ptr<Frame> & frame)
{
  std::cout<<"Cartographer: Add keyframe "<<frame->getFrameId()<<" to GLOBAL map!\n";
  keyframes[frame->getFrameId()] = frame->share_ptr();
  frame->setActive();
  cam_intrinsics[frame->getCamId()] = frame->getCameraIntrinsics();

  // Add new frame to global system
  addFrameToIncSystem();
}


void Cartographer::updateBestMapPointDescriptor(MapLandmark * ml)
{
  // Set best descriptor as the last observation
  // TODO: Decide which is optimal
  MapObservation & mo = (ml->obs_.rbegin())->second;

  feature_extractor_->getDescriptorRaw(mo.frame_ptr->regions_.get(),mo.feat_id,& ml->bestDesc_);
}

void Cartographer::addLandmarksToStructure
(
  Frame * frame,
  std::vector<std::unique_ptr<MapLandmark> > & new_pts,
  const size_t & min_degree_connectivity
)
{
  const size_t frame_id = frame->getFrameId();
  std::cout<<"Cartographer: Add Landmarks to structure: "<<frame->getFrameId()<<": "<<new_pts.size()<<"\n";
  for (std::unique_ptr<MapLandmark> & new_pt : new_pts)
  {
    MapLandmark * m_landmark;
    // Check if point is defined to be put in global map
    // If it has the measurement from this frame we check with adjusted threshold otherwise with a global
    if (new_pt->isValidByConnectivityDegree(new_pt->hasFrameObservation(frame_id) ? min_degree_connectivity : min_obs_per_landmark))
    {
      m_landmark = addLandmarkToGlobalMap(new_pt);
    }
    else
    {
      m_landmark = addLandmarkToLocalMap(new_pt);
      // Mark local landmark as seen in this step
      m_landmark->setObsStep(step_id);
    }

    // Update the best descriptor
    updateBestMapPointDescriptor(m_landmark);

    // Update pointers
    for (auto & obs: m_landmark->obs_)
    {
      MapObservation & mo = obs.second;
      mo.frame_ptr->map_points_[mo.feat_id] = m_landmark;
    }
  }
}

MapLandmark * Cartographer::addLandmarkToLocalMap(std::unique_ptr<MapLandmark> & lm)
{
  MapLandmark * lm_ptr = lm.get();
  tmp_structure[lm_ptr] = std::move(lm);
  return lm_ptr;
}

MapLandmark * Cartographer::addLandmarkToGlobalMap(std::unique_ptr<MapLandmark> & lm)
{
  // Get new ID for global landmark
  const size_t lm_id = getNextFreeLandmarkId();
  MapLandmark * lm_A = lm.get();
  structure[lm_id] = std::move(lm);

  // Update pointers to the landmark
  MapLandmark * lm_tmp = structure[lm_id].get();
  lm_tmp->setActive();

  // Add landmark to system
  addLandmarkToIncSystem();

  return structure[lm_id].get();
}

MapLandmark * Cartographer::addLandmarkFromLocalToGlobalMap(std::unique_ptr<MapLandmark> & lm)
{
  // Get new ID for global landmark
  MapLandmark * g_lm = addLandmarkToGlobalMap(lm);

/*
  for (auto & obs: g_lm->obs_)
  {
   MapObservation & mo = obs.second;
   mo.frame_ptr->map_points_[mo.feat_id] = g_lm;
  }
*/

  // Delete the point from local map
  //tmp_structure.erase(g_lm);

  return g_lm;
}

void Cartographer::addObservationsToLandmarks
(
  Frame * frame,
  const size_t & min_degree_connectivity
)
{
  std::cout<<"Cartographer: Add observations to landmarks! From frame: "<<frame->getFrameId()<<"!\n";
  const size_t & frame_id = frame->getFrameId();
  // Loop through matches and add observation (if valid frame and local point we increase the counter)
  for (size_t feat_id = 0; feat_id < frame->map_points_.size(); ++feat_id)
  {
    MapLandmark * & map_point = frame->map_points_[feat_id];
    if (!map_point)
      continue;
    // Add new observation
    MapObservation m1(feat_id,frame);
    map_point->obs_[frame_id] = m1;

    // Update descriptor to the latest
    updateBestMapPointDescriptor(map_point);

    // If landmark is already in global we add observation to the system
    // If landmark becomes valid with this observation we add it to global map
    if (map_point->isActive())
    {
      addObservationToIncSystem();
    }
    else if (map_point->isValidByConnectivityDegree(min_degree_connectivity))
    {
      // if we have enough valid observations put the landmark to global map
      addLandmarkFromLocalToGlobalMap(tmp_structure[map_point]);
      tmp_structure.erase(map_point);
    }
    else
    {
      // Mark local landmark as seen
      map_point->setObsStep(step_id);
    }
  }
}



bool Cartographer::addStep
(
  std::shared_ptr<Frame> & frame,
  std::vector<std::unique_ptr<MapLandmark> > * vec_new_pts_3D_obs
)
{
  // Increase the step counter
  step_id++;

  const size_t & frame_id = frame->getFrameId();
  size_t min_degree_connectivity = min_obs_per_landmark;

  std::cout<<"Cartographer: Add step "<<frame_id<<"\n";

  // Check if frame is defined in global map
  bool b_frame_def_global = isFrameDefinedInGlobalMap(frame.get());

  // --------------------
  // -- Find minimum degree of connectivity required for the frame to be defined in global map
  // --------------------
  if (!b_frame_def_global)
  {
    min_degree_connectivity = findMinLandmarkDegreeForDefinedInGlobalMap(frame.get(),vec_new_pts_3D_obs);
    if (min_degree_connectivity < 2)
      return false;
  }

  // --------------------
  // -- Add frame to global map
  // --------------------
  addFrameToGlobalMap(frame);
  // --------------------
  // -- Add observations to already existing landmarks
  // --  - If landmark becomes defined in global map we add the whole landmark
  // --------------------
  addObservationsToLandmarks(frame.get(), min_degree_connectivity);

  // --------------------
  // -- Add new points to the global point (if we have sufficient support of valid observations)
  // -- Add new points to local points if we dont have sufficient support
  // --------------------
  if (vec_new_pts_3D_obs)
    addLandmarksToStructure(frame.get(),*vec_new_pts_3D_obs, min_degree_connectivity);

  // Eliminate local points that havent been seen in a long time
  eliminateInactiveLocalLandmarks();

  // Perform global optimization
  optimizeIncSystem();

  return true;
}



// Clears all points and frames in the map (keeps camera intrinsics)
void Cartographer::clearAllMapData()
{
  // Set map to not initialized
  map_initialized = false;

  // Clear everything
  structure.clear();
  tmp_structure.clear();
  // Clear all pointers in keyframes to points in the map
  for (auto kf: keyframes)
  {
    kf.second->clearMapPoints();
  }
  init_map_frames.clear();
  keyframes.clear();
  next_free_landmark_id = 0;

  step_id = 0;

}

void Cartographer::eliminateInactiveLocalLandmarks()
{
  /*
  std::cout<<"Cartographer: Eliminate inactive local landmarks! L_L before: "<<tmp_structure.size()<<"\n";
  for (Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> >::iterator iter_ml = tmp_structure.begin(); iter_ml!=tmp_structure.end();)
  {
    MapLandmark * ml = iter_ml->first;
    if (!ml->isActive())
    {
      if (int(step_id - ml->getLastObsStep()) > max_frames_inactive_local_landmark)
      {
        // To remove local landmark we have to remove any connections to it
        for (auto & obs : ml->obs_)
        {
          obs.second.frame_ptr->clearMapPoint(obs.second.feat_id);
        }
        // remove landmark from local structure
        iter_ml = tmp_structure.erase(iter_ml);

        std::cout<<"Remove local landmark: S: "<<step_id<<" L: "<<ml->getLastObsStep()<<"\n";
      }
      else
      {
        ++iter_ml;
      }
    }
    else
    {
      ++iter_ml;
    }
  }
  std::cout<<" L_L after: "<<tmp_structure.size()<<"\n";*/
}

void Cartographer::getLocalMapPoints(Frame * currentFrame, std::vector<Frame*> & local_frames, std::vector<MapLandmark*> & local_points)
{
  const size_t & frame_cur_id = currentFrame->getFrameId();

  // Get all already reconstructed points from all the frames
  for(Frame * & neigh_frame : local_frames)
  {
    for(MapLandmark * & map_point: neigh_frame->map_points_)
    {
      // If triangulated and not yet accounted for
      if (!map_point || map_point->localMapFrameId_ == frame_cur_id)
        continue;
      local_points.push_back(map_point);
      map_point->localMapFrameId_ = frame_cur_id;
    }
  }
}



/////////////////////////////





















// ------------------------------
// -- Add data to map
// ------------------------------




// ------------------------------
// -- Local map
// ------------------------------


/*
void Cartographer::updateLocalMap(Frame * currentFrame)
{
  local_map_points_.clear();
  // Get neighbour frames of the current frame
  local_map_frames_ = currentFrame->getBestCovisibilityFrames(10);

  // Get all already reconstructed points from all the frames
  for(Frame* neigh_frame : local_map_frames_)
  {
    const size_t frameId = currentFrame->frameId_;
    for(MapLandmark * map_point: neigh_frame->map_points_)
    {
      if (!map_point || map_point->localMapFrameId_ == frameId)
        continue;
      local_map_points_.push_back(map_point);
      map_point->localMapFrameId_ = frameId;
    }
  }
}
*/



// ------------------------------
// -- Camera/Landmark representations
// ------------------------------

Mat34 Cartographer::getCameraProjectionMatrix(Frame * frame, Frame * frame_ref = nullptr)
{
  Mat34 P;
  if (map_camera_type == VSSLAM::MAP_CAMERA_TYPE::GLOBAL)
  {
    return frame->getProjectionMatrix();
  }
  else
  {
    Mat4 T;
    getRelativeCameraTransformation(frame,frame_ref,T);
    return frame->getK()*T.block(0,0,3,4);
  }
}

}
}
