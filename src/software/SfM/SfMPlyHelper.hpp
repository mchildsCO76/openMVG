
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_PLY_HELPER_H
#define OPENMVG_SFM_PLY_HELPER_H

#include "openMVG/numeric/numeric.h"

#include <fstream>
#include <string>
#include <vector>

namespace openMVG{
namespace plyHelper{

/// Export 3D point vector to PLY format
static bool exportToPly(const std::vector<Vec3> & vec_points,
  const std::string & sFileName)
{
  std::ofstream outfile(sFileName.c_str());
  if (!outfile.is_open())
    return false;

  outfile << "ply"
    << std::endl << "format ascii 1.0"
    << std::endl << "element vertex " << vec_points.size()
    << std::endl << "property float x"
    << std::endl << "property float y"
    << std::endl << "property float z"
    << std::endl << "property uchar red"
    << std::endl << "property uchar green"
    << std::endl << "property uchar blue"
    << std::endl << "end_header" << std::endl;

  for (size_t i=0; i < vec_points.size(); ++i)
  {
    outfile
      << vec_points[i](0) << ' '
      << vec_points[i](1) << ' '
      << vec_points[i](2) << ' '
      << "255 255 255" << "\n";
  }
  const bool bOk = outfile.good();
  outfile.close();
  return bOk;
}

/// Export 3D point vector and camera position to PLY format
static bool exportToPly(const std::vector<Vec3> & vec_points,
  const std::vector<Vec3> & vec_camPos,
  const std::string & sFileName,
  const std::vector<Vec3> * vec_coloredPoints = NULL,
  const std::vector<double> * vec_structureUncertainty = NULL)
{
  std::ofstream outfile(sFileName.c_str());
  if (!outfile.is_open())
    return false;


  outfile << "ply"
    << '\n' << "format ascii 1.0"
    << '\n' << "element vertex " << vec_points.size()+vec_camPos.size()
    << '\n' << "property float x"
    << '\n' << "property float y"
    << '\n' << "property float z"
    << '\n' << "property uchar red"
    << '\n' << "property uchar green"
    << '\n' << "property uchar blue";

  if(vec_structureUncertainty!=NULL){
    outfile << '\n' << "property float quality";
  }

  outfile << '\n' << "end_header" << std::endl;

  for (size_t i=0; i < vec_points.size(); ++i)  {
    outfile
      << vec_points[i](0) << ' '
      << vec_points[i](1) << ' '
      << vec_points[i](2) << ' ';

    if (vec_coloredPoints == NULL){
      outfile << "255 255 255";
    }else{
      outfile
          << (*vec_coloredPoints)[i](0) << ' '
          << (*vec_coloredPoints)[i](1) << ' '
          << (*vec_coloredPoints)[i](2);
    }

    if(vec_structureUncertainty == NULL){
      outfile<<'\n';
    }
    else{
      outfile << ' ' << (*vec_structureUncertainty)[i]<<'\n';
    }

  }

  for (size_t i=0; i < vec_camPos.size(); ++i)  {
    outfile
      << vec_camPos[i](0) << ' '
      << vec_camPos[i](1) << ' '
      << vec_camPos[i](2) << ' '
      << "0 255 0 0.0 " << "\n";
  }
  outfile.flush();
  const bool bOk = outfile.good();
  outfile.close();
  return bOk;
}


} // namespace plyHelper
} // namespace openMVG

#endif // OPENMVG_SFM_PLY_HELPER_H

