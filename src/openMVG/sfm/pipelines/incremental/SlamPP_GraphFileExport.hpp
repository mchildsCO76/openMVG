
#ifndef OPENMVG_SLAMPP_GRAPH_FILE_EXPORT
#define OPENMVG_SLAMPP_GRAPH_FILE_EXPORT

#include <cstdlib>

#include <openMVG/types.hpp>
#include "openMVG/sfm/sfm.hpp"

#include "openMVG/sfm/pipelines/incremental/incremental_SfM.hpp"

#include "lemon/smart_graph.h"

namespace openMVG {
namespace sfm {

class IncrementalSfMReconstructionEngine;

struct SlamPP_Data
{
  IndexT next_free_id_slamPP = 0;

  /// SlamPP index mapping
  Hash_Map<IndexT,IndexT> camera_ids_omvg_slamPP;
  Hash_Map<IndexT,IndexT> track_ids_omvg_slamPP;
  Hash_Map<IndexT,IndexT> camera_ids_slamPP_omvg;
  Hash_Map<IndexT,IndexT> track_ids_slamPP_omvg;

  // Owner of tracks
  Hash_Map<IndexT,IndexT> owner_track_cam_id;
  // Parents graph
  lemon::SmartGraph g_parents_cam;
  std::unique_ptr<lemon::SmartGraph::NodeMap<IndexT> > g_parents_nodeMap;
  Hash_Map<IndexT, lemon::SmartGraph::Node > g_parents_indexMap;
  //std::unique_ptr<lemon::SmartGraph::NodeMap<IndexT> > map_g_node_camId_omvg;
  //Hash_Map<IndexT, lemon::SmartGraph::Node> map_g_camId_omvg_node;
  std::unique_ptr<lemon::Bfs<lemon::SmartGraph> > bfs;
  

  /// Export format
  int iOutputVertexType = 0; //0 - SE3; 1 - Sim3
  bool iOutputLandmarkType = 0; // 0 - Eucliean (world); 1 - inverse depth (reference)

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

  void setGraphFileOutputFile(std::string &filename)
  {
    slamPP_dataset_filename = filename;
  }

  void initCamParentsGraph()
  {
    g_parents_nodeMap.reset(new lemon::SmartGraph::NodeMap<IndexT>(g_parents_cam) );
  }

  void addCamWithParentToGraph(IndexT camId_omvg, IndexT cam_parentId_omvg)
  {
    // Add new camera to the graph
    std::cout<<"Add camera "<<camId_omvg<<" to graph\n";
    lemon::SmartGraph::Node n = g_parents_cam.addNode();
    (*g_parents_nodeMap)[n] = camId_omvg;
    g_parents_indexMap[camId_omvg] = n;

    // Create edge between the camera and its parent if there is a parent
    if (cam_parentId_omvg != std::numeric_limits<IndexT>::max())
    {
      std::cout<<"Add LINK "<<camId_omvg<< " :: "<<cam_parentId_omvg<<" to graph\n";
      lemon::SmartGraph::Node sourceNode = g_parents_indexMap[camId_omvg];
      lemon::SmartGraph::Node targerNode = g_parents_indexMap[cam_parentId_omvg];
      g_parents_cam.addEdge(sourceNode,targerNode);
    }
  }

  
  void initBreadthSearchFirstGraph()
  {
    bfs.reset(new lemon::Bfs<lemon::SmartGraph>(g_parents_cam));
  }
  
  void findPathBetween(IndexT camA, IndexT camB, std::vector<IndexT> &cam_predecessor)
  {
    lemon::SmartGraph::Node startN = g_parents_indexMap[camA];
    lemon::SmartGraph::Node endN = g_parents_indexMap[camB];
    bfs->run(startN,endN);
    for (lemon::SmartGraph::Node v = endN; v != startN; v = bfs->predNode(v))
    {
      if (v != lemon::INVALID && bfs->reached(v)) //special LEMON node constant
      {
         cam_predecessor.push_back((*g_parents_nodeMap)[v]);
      }
    }
  }

  SlamPP_Data()
  {
   initCamParentsGraph();
  }

};


} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SLAMPP_GRAPH_FILE_EXPORT
