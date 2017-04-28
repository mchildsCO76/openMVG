// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>
#include <iostream>


#include "openMVG/image/image_io.hpp"
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <openMVG/vsslam/system/Camera.hpp>
#include <openMVG/vsslam/system/slam_system.hpp>
#include <openMVG/vsslam/display/vsslam_display.hpp>

#include <openMVG/vsslam/tracking/Abstract_Tracker.hpp>
#include <openMVG/vsslam/tracking/Tracker_Features.hpp>
#include <openMVG/vsslam/features/Abstract_Feature_Extractor.hpp>
#include <openMVG/vsslam/features/Feat_Extractor_SIFT.hpp>
#include <openMVG/vsslam/features/Abstract_Feature_Matcher.hpp>
#include <openMVG/vsslam/features/Feat_Matcher_CascadeHashing.hpp>
#include <openMVG/vsslam/features/Feat_Matcher_Regions.hpp>

#include <openMVG/vsslam/mapping/Cartographer.hpp>


#include "software/VSSLAM/CGlWindow.hpp"

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::image;
using namespace openMVG::vsslam;

VSSLAM_Display display_data;
VSSLAM_Time_Stats time_data;

int n_dummy_param = 0;

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "VISUAL S-SLAM -- Tracking demo --" << std::endl;

  CmdLine cmd;

  std::string sImaDirectory = "";
  std::string sImaMask = "";
  unsigned int uTracker = 0;

  // Camera data
  std::string sKmatrix;
  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('m', sImaMask, "imamask") );
  cmd.add( make_option('t', uTracker, "tracker") );
  cmd.add( make_option('k', sKmatrix, "intrinsics") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--imadir path] \n"
    << "[-m|--mask image] \n"
    << "[-t|--tracker Used tracking interface] \n"
    << "\t 0 (default) description based Tracking -> Fast detector + Dipole descriptor\n"
#if defined HAVE_OPENCV
    << "\t 1 image based Tracking -> use OpenCV Pyramidal KLT Tracking\n"
#endif
    << "[-k|--intrinsics] Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\"\n"
    << "[-c|--camera_model] Camera model type:\n"
    << "\t 1: Pinhole\n"
    << "\t 2: Pinhole radial 1\n"
    << "\t 3: Pinhole radial 3 (default)\n"
    << "\t 4: Pinhole brown 2\n"
    << "\t 5: Pinhole with a simple Fish-eye distortion\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

   std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImaDirectory << std::endl
            << "--imageMask " << sImaMask << std::endl
            << "--intrinsics " << sKmatrix << std::endl
            << "--camera_model " << i_User_camera_model << std::endl;

  if (sImaDirectory.empty() || !stlplus::is_folder(sImaDirectory))
  {
    std::cerr << "\nIt is an invalid input directory" << std::endl;
    return EXIT_FAILURE;
  }


  // ----------------------------------
  // Image management
  // ----------------------------------
  image::Image<unsigned char> currentImage;
  // Load images from folder
  std::vector<std::string> vec_image = stlplus::folder_files(sImaDirectory);
  // clean invalid image file
  {
    std::vector<std::string> vec_image_;
    for (size_t i = 0; i < vec_image.size(); ++i)
    {
      if(vec_image[i].find("mask.png") != std::string::npos
         || vec_image[i].find("_mask.png") != std::string::npos)
      {
        std::cout
            << vec_image[i] << " is a mask image" << "\n";
        continue;
      }

      if (openMVG::image::GetFormat(vec_image[i].c_str()) != openMVG::image::Unknown)
        vec_image_.push_back(vec_image[i]);
    }
    vec_image_.swap(vec_image);
  }
  std::sort(vec_image.begin(), vec_image.end());

  // ----------------------------------
  // SLAM system initialization
  // ----------------------------------
  std::cout<<"\nVSSLAM: [Start]\n";
  std::shared_ptr<VSSLAM_Parameters> params_system = std::make_shared<VSSLAM_Parameters>();

  SLAM_System slam_system(params_system);

  // Create cartographer
  MAP_FRAME_TYPE map_frame_type = MAP_FRAME_TYPE::GLOBAL;
  MAP_LANDMARK_TYPE map_landmark_type = MAP_LANDMARK_TYPE::GLOBAL_EUCLIDEAN;

  MAP_OPTIMIZATION_TYPE global_BA_type = MAP_OPTIMIZATION_TYPE::SLAMPP;
  MAP_OPTIMIZATION_TYPE local_BA_type = MAP_OPTIMIZATION_TYPE::CERES;

  if (!slam_system.createCartographer(map_frame_type,map_landmark_type,global_BA_type,local_BA_type ))
  {
    std::cerr << "Cannot instantiate the cartographer" << std::endl;
    return EXIT_FAILURE;
  }


  {
    std::unique_ptr<Abstract_Tracker> ptr_tracker;
    std::unique_ptr<Abstract_Feature_Extractor> ptr_feat_extractor;
    std::unique_ptr<Abstract_Feature_Matcher> ptr_feat_matcher;

    switch (uTracker)
    {
      case 0:
        ptr_feat_extractor.reset(new Feat_Extractor_SIFT(params_system, features::HIGH_PRESET));
        ptr_feat_matcher.reset(new Feat_Matcher_CascadeHashing(params_system, ptr_feat_extractor.get()));
        display_data.b_enable_display = 0;

        break;
      case 1:
        ptr_feat_extractor.reset(new Feat_Extractor_SIFT(params_system, features::HIGH_PRESET));
        ptr_feat_matcher.reset(new Feat_Matcher_CascadeHashing(params_system, ptr_feat_extractor.get()));
        display_data.b_enable_display = 1;
        //ptr_feat_matcher.reset(new Feat_Matcher_Regions(params_system, ptr_feat_extractor.get()));

        break;
      default:
        std::cerr << "Unknow tracking method" << std::endl;
        return EXIT_FAILURE;
    }

    ptr_tracker.reset(new Tracker_Features(params_system));
    if (!ptr_tracker)
    {
      std::cerr << "Cannot instantiate the tracking interface" << std::endl;
      return EXIT_FAILURE;
    }

    slam_system.setFeatureExtractor(ptr_feat_extractor);
    slam_system.setFeatureMatcher(ptr_feat_matcher);
    slam_system.setTracker(ptr_tracker);


  }



  // ----------------------------------
  // Create camera
  // ----------------------------------
  IndexT id_cam_0;
  {
    // Check and set camera parameters
    CameraParameters params_cam_0;
    params_cam_0.camera_model = EINTRINSIC(i_User_camera_model);

    if (sKmatrix.empty() ||
      !CameraParameters::checkIntrinsicStringValidity(sKmatrix, params_cam_0.focal, params_cam_0.ppx, params_cam_0.ppy) )
    {
      std::cerr << "\nInvalid K matrix input" << std::endl;
      return EXIT_FAILURE;
    }
    if (!params_cam_0.readImageSettings(stlplus::create_filespec( sImaDirectory, vec_image[0] )))
    {
      std::cerr << "\nError reading image header file" << std::endl;
      return EXIT_FAILURE;
    }

    // Create camera model
    if (params_cam_0.isValid())
    {
      params_cam_0.b_calibrated = true;
      id_cam_0 = slam_system.createCamera(params_cam_0);
    }
    else
    {
      std::cerr << "Error: invalid camera parameters"<< std::endl;
      return EXIT_FAILURE;
    }
  }
  // ----------------------------------
  // Load mask images
  // ----------------------------------
  {
    image::Image<unsigned char> mask_cam_0;
    if (!sImaMask.empty() && stlplus::file_exists(sImaMask))
    {
      if (image::GetFormat(sImaMask.c_str()) == image::Unknown)
        std::cout << "\nMask image path is invalid! Not using mask image!" << std::endl;
      else
      {
        if (!(image::ReadImage( sImaMask.c_str(), &mask_cam_0)))
        {
          std::cout << "\nMask image is invalid! Not using mask image!" << std::endl;
        }
        else
        {
          std::cout << "\nUsing mask image: "<<sImaMask<<"!\n" << std::endl;
          slam_system.addMaskImageToCamera(id_cam_0,mask_cam_0);
        }
      }
    }
  }

  if(!slam_system.isReady())
  {
    std::cerr << "VSSLAM: SLAM not correctly initialized\n";
    return EXIT_FAILURE;
  }

  // ----------------------------------
  // Graphics
  // ----------------------------------
  if ( !glfwInit() )
  {
    return EXIT_FAILURE;
  }

  CGlWindow window;
  GLuint text2D;

  // ----------------------------------
  // Frame-by-Frame processing
  // ----------------------------------
  IndexT id_frame = 0;
  IndexT id_cam = 0;
  double timestamp_frame = 0;

  // Traverse image folder
  for (std::vector<std::string>::const_iterator iterFile = vec_image.begin();
      iterFile != vec_image.end(); ++iterFile, ++id_frame, ++timestamp_frame)
  {
    const std::string sImageFilename = stlplus::create_filespec( sImaDirectory, *iterFile );
    if (openMVG::image::ReadImage( sImageFilename.c_str(), &currentImage))
    {
      slam_system.nextFrame(currentImage, id_frame, id_cam,timestamp_frame);
std::cout<<"ID FRAME: "<<id_frame<<"\n";
      if (window._height < 0)
      {
        // no window created yet, initialize it with the first frame

        const double aspect_ratio = currentImage.Width()/(double)currentImage.Height();
        window.Init(1280, 1280/aspect_ratio, "VisualOdometry--TrackingViewer");
        glGenTextures(1,&text2D);             //allocate the memory for texture
        glBindTexture(GL_TEXTURE_2D, text2D); //Binding the texture
        glEnable(GL_TEXTURE_2D);              //Enable texture
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE,
          currentImage.Width(), currentImage.Height(), 0,
          GL_LUMINANCE, GL_UNSIGNED_BYTE, currentImage.data());
      }

      if (display_data.b_enable_display)
      {
        // Display steps
        display_data.displaySteps(window,text2D, currentImage,slam_system.getCurrentFramePtr(),2);

        display_data.resetSteps();
      }
      else
      {
        display_data.displayImage(window,text2D, currentImage);
        display_data.displayDetectedFeatures(slam_system.getCurrentFramePtr());
        display_data.displayHistoryTracks(slam_system.getCurrentFramePtr());
      }

      glFlush();
      window.Swap(); // Swap openGL buffer
      std::cout<<"Press ENTER to continue....."<<std::endl<<std::endl;
      //std::cin.ignore(1);
    }
  }





  return 0;
}
