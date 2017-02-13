
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

  double img_border[4]; // borders of image (min_x,min_y,max_x,max_y)

  void computeImageBorders()
  {
    if (bCalibrated && cam_intrinsic_->have_disto())
    {
      // If camera is calibrated and original intrinsic have distortion we compute the new borders by undistorting control points

      // Control points:
      // 0 -- 4 -- 1
      // |         |
      // 7         5
      // |         |
      // 3 -- 6 -- 2
      const size_t w = cam_intrinsic_->w();
      const size_t h = cam_intrinsic_->h();
      Mat2X pt2D_cp(2, 8);
      pt2D_cp << 0,w,w,0,w/2,w,w/2,0,     0, 0, h, h, 0, h/2, h, h/2;
      std::cout<<"CPP: "<<pt2D_cp<<"\n";
      cam_intrinsic_undist->get_ud_pixel(pt2D_cp);
      std::cout<<"CP: "<<pt2D_cp<<"\n";


    }
    else
    {
      // Use image corners as limits as we will always do all the processing on distorted coordinates
      img_border[0] = 0;
      img_border[1] = 0;
      img_border[2] = cam_intrinsic_ptr->w();
      img_border[3] = cam_intrinsic_ptr->h();
    }
  }

  bool isPointInImageBorders(const Vec2 & pt) const
  {
    if (pt(0) < img_border[0] || pt(0) > img_border[2] || pt(1) < img_border[1] || pt(1) > img_border[3])
      return false;
    return true;
  }
};

}
}
