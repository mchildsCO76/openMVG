
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
  std::vector<std::unique_ptr<MapLandmark> > * vec_new_pts_3D
)
{
  step_id++;

  const IndexT & frame_id = frame->getFrameId();

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
    if (vec_new_pts_3D)
      addLandmarksToStructure(frame.get(),*vec_new_pts_3D, min_obs_per_landmark+1);

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
      for (auto & frame_init_map : init_map_frames)
      {
        addFrameToGlobalMap(frame_init_map);
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
  std::vector<std::unique_ptr<MapLandmark> > * vec_new_pts_3D
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
  if (vec_new_pts_3D)
  {
    const size_t frame_id = frame->getFrameId();
    for (std::unique_ptr<MapLandmark> & ml : *vec_new_pts_3D)
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

  feature_extractor_->getDescriptorRaw(mo.frame_ptr->getRegions(),mo.feat_id,& ml->bestDesc_);
}

void Cartographer::addLandmarksToStructure
(
  Frame * frame,
  std::vector<std::unique_ptr<MapLandmark> > & new_3D_pts,
  const size_t & min_degree_connectivity
)
{
  const IndexT frame_id = frame->getFrameId();
  std::cout<<"Cartographer: Add Landmarks to structure: "<<frame->getFrameId()<<": "<<new_3D_pts.size()<<"\n";
  for (std::unique_ptr<MapLandmark> & pt_3D_new : new_3D_pts)
  {
    MapLandmark * m_landmark;
    // Check if point is defined to be put in global map
    // If it has the measurement from this frame we check with adjusted threshold otherwise with a global
    if (map_initialized && pt_3D_new->isValidByConnectivityDegree(pt_3D_new->hasFrameObservation(frame_id) ? min_degree_connectivity : min_obs_per_landmark))
    {
      m_landmark = addLandmarkToGlobalMap(pt_3D_new);
    }
    else
    {
      m_landmark = addLandmarkToLocalMap(pt_3D_new);
      // Mark local landmark as seen in this step
      m_landmark->setObsStep(step_id);
    }

    // Update the best descriptor
    updateBestMapPointDescriptor(m_landmark);


    // Create connections between frame and map_point
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
  lm_ptr->id_ = local_p_id++;
  // Update count of points for each frame statistics
  for (auto & obs: lm_ptr->obs_)
  {
    increaseLocalFrameCount(obs.second.frame_ptr);
  }

  return lm_ptr;
}

MapLandmark * Cartographer::addLandmarkToGlobalMap(std::unique_ptr<MapLandmark> & lm)
{
  // Get new ID for global landmark
  const IndexT lm_id = getNextFreeLandmarkId();
  MapLandmark * lm_A = lm.get();
  structure[lm_id] = std::move(lm);

  // Update pointers to the landmark
  MapLandmark * lm_tmp = structure[lm_id].get();
  lm_tmp->setActive();

  // Add landmark to system
  addLandmarkToIncSystem();

  return lm_tmp;
}

MapLandmark * Cartographer::addLandmarkFromLocalToGlobalMap(std::unique_ptr<MapLandmark> & lm)
{
  // Update the statistics in local frames
  // Check if frame is already in local frames
  for (auto & obs: lm->obs_)
  {
    decreaseLocalFrameCount(obs.second.frame_ptr);
  }


  // Get new ID for global landmark
  MapLandmark * g_lm = addLandmarkToGlobalMap(lm);

  return g_lm;
}

void Cartographer::addObservationsToLandmarks
(
  Frame * frame,
  const size_t & min_degree_connectivity
)
{
  std::cout<<"Cartographer: Add observations to landmarks! From frame: "<<frame->getFrameId()<<"!\n";
  const IndexT & frame_id = frame->getFrameId();
  // Loop through matches and add observation (if valid frame and local point we increase the counter)
  for (IndexT feat_id = 0; feat_id < frame->map_points_.size(); ++feat_id)
  {
    MapLandmark * & map_point = frame->map_points_[feat_id];
    if (!map_point)
      continue;

    // If landmark is already in global we add observation to the system
    // If landmark becomes valid with this observation we add it to global map
    if (map_point->isActive())
    {
      // Add observation to the landmark
      map_point->addObservation(frame,feat_id);
      addObservationToIncSystem();
    }
    else if (map_point->isValidByConnectivityDegree(min_degree_connectivity))
    {
      // if we have enough valid observations put the landmark to global map
      map_point = addLandmarkFromLocalToGlobalMap(tmp_structure[map_point]);

      // Add observation to the landmark
      map_point->addObservation(frame,feat_id);
      // Erase the point from local map
      tmp_structure.erase(map_point);

      addLandmarkToIncSystem();

    }
    else
    {
      // Add observation to the landmark
      map_point->addObservation(frame,feat_id);
      // Mark local landmark as seen
      map_point->setObsStep(step_id);

      // Update the statistics in local frames
      increaseLocalFrameCount(frame);
    }

    // Update descriptor to the latest
    updateBestMapPointDescriptor(map_point);

  }

}

void Cartographer::verifyLocalLandmarks(Frame * frame)
{
  // Go through local frames and check the ones which are not fixed
  std::vector<size_t> vec_tmp_structure_outliers;
  #ifdef OPENMVG_USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for(size_t t_i = 0; t_i < tmp_structure.size(); ++t_i)
  {
    Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> >::iterator it_tmp_structure = tmp_structure.begin();
    std::advance(it_tmp_structure, t_i);
    MapLandmark * ml = it_tmp_structure->first;

    LandmarkObservations & obs = ml->obs_;

    Vec3 pt_3D_frame_i;
    for(LandmarkObservations::iterator iter_mo = obs.begin(); iter_mo != obs.end();)
    {
      MapObservation & m_o =  iter_mo->second;
      Frame * & frame_i = m_o.frame_ptr;
      // Put 3D point into coordinate system of frame
      getRelativePointPosition(ml->X_,ml->ref_frame_,pt_3D_frame_i,frame_i);

      if (frame_i->getSquaredReprojectionError(pt_3D_frame_i,m_o.feat_id) > frame_i->AC_reprojection_thresh_)
      {
        std::cout<<" Remove local point - measurement: "<<frame_i->getFrameId()<<":: "<<t_i<<" :: "<<ml->id_<<"\n";
        // Remove 3D point from frame
       frame_i->clearMapPoint(m_o.feat_id);
        // Update local frame count
        decreaseLocalFrameCount(frame_i);
        iter_mo = obs.erase(iter_mo);
      }
      else
      {
        ++iter_mo;
      }
    }

    if (obs.size()<2)
    {
      if (obs.size() == 1)
      {
        // Delete the 3D point from remaining frame
        MapObservation & m_o =  obs.begin()->second;
        std::cout<<" Remove local point - remaining measurement: "<<m_o.frame_ptr->getFrameId()<<" :: "<<ml->id_<<"\n";
        m_o.frame_ptr->clearMapPoint(m_o.feat_id);
        decreaseLocalFrameCount(m_o.frame_ptr);
      }
      // Delete point from the vector of new points
      #pragma omp critical
      {
        vec_tmp_structure_outliers.push_back(t_i);
      }
    }

  }

  // sort indexes of outliers
  std::sort(vec_tmp_structure_outliers.begin(),vec_tmp_structure_outliers.end());

  // Remove any triangulated landmarks that dont have enough measurements
  for (size_t o_i = 0; o_i < vec_tmp_structure_outliers.size(); ++o_i)
  {
    Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> >::iterator it_outlier = tmp_structure.begin();
    // Reduce the number by the number of elements already deleted
    std::advance(it_outlier,vec_tmp_structure_outliers[o_i]-o_i);

    std::cout<<" Remove local point - outlier :: "<<it_outlier->first->obs_.size()<<": "<<o_i<<" :: "<<vec_tmp_structure_outliers[o_i]-o_i<<" :: "<<it_outlier->first->id_<<"\n";

    // Remove possible connection with current local frame
    for (IndexT feat_id = 0; feat_id < frame->map_points_.size(); ++feat_id)
    {
      MapLandmark * & ml = frame->map_points_[feat_id];
      if (ml && ml == it_outlier->first)
      {
        std::cout<<"FF: "<<feat_id<<"\n";
        frame->clearMapPoint(feat_id);
        break;
      }
    }

    // Delete 3D point (all other references have been deleted before)
    std::cout<<"TTB: "<<tmp_structure.size()<<"\n";
    tmp_structure.erase(it_outlier);
    std::cout<<"TTA: "<<tmp_structure.size()<<"\n";
  }
}

void Cartographer::decreaseLocalFrameCount(Frame * frame)
{
  Hash_Map<Frame*, size_t>::iterator it_tmp_frames = tmp_frames.find(frame);
  if (it_tmp_frames->second>1)
    it_tmp_frames->second--;
  else
    tmp_frames.erase(it_tmp_frames);
}
void Cartographer::increaseLocalFrameCount(Frame * frame)
{
  Hash_Map<Frame*, size_t>::iterator it_tmp_frames = tmp_frames.find(frame);
  if (it_tmp_frames==tmp_frames.end())
    tmp_frames[frame] = 1;
  else
    it_tmp_frames->second++;
}

bool Cartographer::addStep
(
  std::shared_ptr<Frame> & frame,
  std::vector<std::unique_ptr<MapLandmark> > * vec_new_pts_3D
)
{
  // Increase the step counter
  step_id++;

  const IndexT & frame_id = frame->getFrameId();
  size_t min_degree_connectivity = min_obs_per_landmark;

  std::cout<<"Cartographer: Add step "<<frame_id<<"\n";

  // Check if frame is defined in global map
  bool b_frame_def_global = isFrameDefinedInGlobalMap(frame.get());

  // --------------------
  // -- Find minimum degree of connectivity required for the frame to be defined in global map
  // --------------------
  if (!b_frame_def_global)
  {
    min_degree_connectivity = findMinLandmarkDegreeForDefinedInGlobalMap(frame.get(),vec_new_pts_3D);
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
  if (vec_new_pts_3D)
    addLandmarksToStructure(frame.get(),*vec_new_pts_3D, min_degree_connectivity);

  // Eliminate local points that havent been seen in a long time
  double start_time = omp_get_wtime();
  eliminateInactiveLocalLandmarks();
  std::cout<<"Eliminate: "<<omp_get_wtime()-start_time<<"\n";
  // Perform global optimization
  optimizeIncSystem();
std::cout<<"MAP STATUS: P: "<<structure.size()<<" L: "<<tmp_structure.size()<<" K: "<<keyframes.size()<<"\n";

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
  std::cout<<"Cartographer: Eliminate inactive local landmarks! Before: "<<tmp_structure.size();
  for (Hash_Map<MapLandmark*, std::unique_ptr<MapLandmark> >::iterator iter_ml = tmp_structure.begin(); iter_ml!=tmp_structure.end();)
  {
    MapLandmark * ml = iter_ml->first;

    if (int(step_id - ml->getLastObsStep()) > max_frames_inactive_local_landmark)
    {
      // To remove local landmark we have to remove any connections to it
      for (auto & obs : ml->obs_)
      {
        obs.second.frame_ptr->clearMapPoint(obs.second.feat_id);

        // Update the statistics in local frames
        decreaseLocalFrameCount(obs.second.frame_ptr);
      }
      // remove landmark from local structure
      iter_ml = tmp_structure.erase(iter_ml);

    }
    else
    {
      ++iter_ml;
    }
  }
  std::cout<<" After: "<<tmp_structure.size()<<"\n";
}

void Cartographer::getLocalMapPoints(Frame * frame_current, std::vector<Frame*> & local_frames, std::vector<MapLandmark*> & local_points)
{
  const IndexT & frame_cur_id = frame_current->getFrameId();

  // Get all already reconstructed points from all the frames
  for(Frame * & frame_neigh : local_frames)
  {
    for(MapLandmark * & map_point: frame_neigh->map_points_)
    {
      // If triangulated and not yet accounted for
      if (!map_point || map_point->last_local_map_frame_id_ == frame_cur_id)
        continue;
      local_points.push_back(map_point);
      map_point->last_local_map_frame_id_ = frame_cur_id;
    }
  }
}

bool Cartographer::exportSceneToPly(const std::string & filename,
    sfm::ESfM_Data flags_part)
{
  const bool b_structure = (flags_part & sfm::STRUCTURE) == sfm::STRUCTURE;
  const bool b_local_structure = (flags_part & sfm::CONTROL_POINTS) == sfm::CONTROL_POINTS;
  const bool b_extrinsics = (flags_part & sfm::EXTRINSICS) == sfm::EXTRINSICS;

  // Create the stream and check its status
  std::ofstream stream(filename.c_str());
  if (!stream.is_open())
    return false;

  bool bOk = false;
  // Count how many views having valid poses:
  IndexT keyframes_count = 0;
  if (b_extrinsics)
  {
    keyframes_count = keyframes.size();
  }

  stream << std::fixed << std::setprecision (std::numeric_limits<double>::digits10 + 1);

  stream << "ply"
    << '\n' << "format ascii 1.0"
    << '\n' << "element vertex "
      // Vertex count: (#landmark + #GCP + #view_with_valid_pose)
      << (  (b_structure ? structure.size() : 0)
          + (b_local_structure ? tmp_structure.size() : 0)
          + keyframes_count)
    << '\n' << "property double x"
    << '\n' << "property double y"
    << '\n' << "property double z"
    << '\n' << "property uchar red"
    << '\n' << "property uchar green"
    << '\n' << "property uchar blue"
    << '\n' << "end_header" << std::endl;


  if (b_extrinsics)
  {
    for (auto & frame_it : keyframes)
    {
      // Export pose as Green points
      Vec3 center = frame_it.second->getCameraCenter();
        stream
          << center(0) << ' '
          << center(1) << ' '
          << center(2) << ' '
          << "0 255 0\n";
    }
  }

  if (b_structure)
  {
    // Export structure points as White points
    for ( const auto & iterMapLandmarks : structure )
    {
      MapLandmark * map_point = iterMapLandmarks.second.get();
      Vec3 pt_3D_w;
      getRelativePointPosition(map_point->X_,map_point->ref_frame_,pt_3D_w,nullptr);
      stream
        << pt_3D_w(0) << ' '
        << pt_3D_w(1) << ' '
        << pt_3D_w(2) << ' '
        << "255 255 255\n";
    }
  }


  if (b_local_structure)
  {
    // Export local structure points as red points
    for ( const auto & iterMapLandmarks : tmp_structure )
    {
      MapLandmark * map_point = iterMapLandmarks.second.get();
      Vec3 pt_3D_w;
      getRelativePointPosition(map_point->X_,map_point->ref_frame_,pt_3D_w,nullptr);
      stream
        << pt_3D_w(0) << ' '
        << pt_3D_w(1) << ' '
        << pt_3D_w(2) << ' '
        << "255 0 0\n";
    }
  }

  stream.flush();
  bOk = stream.good();
  stream.close();

  return bOk;

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
