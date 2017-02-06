
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <openMVG/cameras/cameras.hpp>
#include <iostream>
#include <memory>
#include "openMVG/image/image.hpp"


namespace openMVG  {
namespace VSSLAM  {

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;


struct CameraParams
{
  EINTRINSIC camera_model;
  double focal = -1;
  double ppx = -1;
  double ppy = -1;
  size_t img_width = 0;
  size_t img_height = 0;
  std::vector<double> dist;
  bool bCalibrated = true;

  bool readImageSettings(const std::string & sImagePath)
  {
    ImageHeader imgHeader;

    if (!openMVG::image::ReadImageHeader(sImagePath.c_str(), &imgHeader))
    {
      std::cerr << "\nError reading image header file" << std::endl;
      return false;
    }
    img_width = imgHeader.width;
    img_height = imgHeader.height;
    // If principal point is not set we set the center of image
    if (ppx == -1)
      ppx = img_width/2;
    if (ppy == -1)
      ppy = img_height/2;

    return true;
  }

  bool checkValidParams()
  {
    if (focal !=-1 && ppx != -1 && ppy != -1 && img_width != 0 && img_height != 0)
    {
      return true;
    }
    return false;
  }

};

/// Camera
class Camera : public std::enable_shared_from_this<Camera>
{
public:
  size_t cam_id;
  std::shared_ptr<IntrinsicBase> cam_intrinsic_;
  std::shared_ptr<IntrinsicBase> cam_intrinsic_undist;
  IntrinsicBase * cam_intrinsic_ptr = nullptr; // Pointitng to the camera that should be used in processing (distorted/undistorted)
  bool bCalibrated = false;
};

}
}
