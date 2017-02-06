
#pragma once

// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <openMVG/sfm/sfm_data_io.hpp>

#include <fstream>
#include <iomanip>
using namespace openMVG;


namespace openMVG {
namespace VSSLAM {

/// Save the structure and camera positions of a SfM_Data container as 3D points in a PLY ASCII file.
inline bool Save_PLY(
  const VSSLAM_Data & slam_data,
  const std::string & filename,
  sfm::ESfM_Data flags_part)
{
  const bool b_structure = (flags_part & sfm::ESfM_Data::STRUCTURE) == sfm::ESfM_Data::STRUCTURE;
  const bool b_extrinsics = (flags_part & sfm::ESfM_Data::EXTRINSICS) == sfm::ESfM_Data::EXTRINSICS;

  if (!(b_structure || b_extrinsics))
    return false; // No 3D points to display, so it would produce an empty PLY file

  // Create the stream and check its status
  std::ofstream stream(filename.c_str());
  if (!stream.is_open())
    return false;

  bool bOk = false;
  {
    // Count how many views having valid poses:
    IndexT keyframe_count = b_extrinsics ? slam_data.keyframes.size(): 0;

    stream << std::fixed << std::setprecision (std::numeric_limits<double>::digits10 + 1);

    stream << "ply"
      << '\n' << "format ascii 1.0"
      << '\n' << "element vertex "
        // Vertex count: (#landmark + #GCP + #view_with_valid_pose)
        << (  (b_structure ? slam_data.structure.size() : 0)
            + keyframe_count)
      << '\n' << "property double x"
      << '\n' << "property double y"
      << '\n' << "property double z"
      << '\n' << "property uchar red"
      << '\n' << "property uchar green"
      << '\n' << "property uchar blue"
      << '\n' << "end_header" << std::endl;

      if (b_extrinsics)
      {
        for (const auto & keyframe : slam_data.keyframes)
        {
          // Export pose as Green points

            stream
              << keyframe.second->pose_.center()(0) << ' '
              << keyframe.second->pose_.center()(1) << ' '
              << keyframe.second->pose_.center()(2) << ' '
              << "0 255 0\n";
        }
      }

      if (b_structure)
      {
        // Export structure points as White points
        for ( const auto & iterLandmarks : slam_data.structure )
        {
          stream
            << iterLandmarks.second.X_(0) << ' '
            << iterLandmarks.second.X_(1) << ' '
            << iterLandmarks.second.X_(2) << ' '
            << "255 255 255\n";
        }
      }

      stream.flush();
      bOk = stream.good();
      stream.close();
  }
  return bOk;
}

}
}
