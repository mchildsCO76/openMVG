
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/features/feature.hpp"
#include <software/VO/CGlWindow.hpp>

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include <stdlib.h>
#include <iostream>

// Tracker
#include <software/VSSLAM/slam/Abstract_Tracker.hpp>
#include <software/VSSLAM/slam/Abstract_FeatureExtractor.hpp>

using namespace openMVG;

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "VISUAL SSLAM -- Tracking demo --" << std::endl;

  CmdLine cmd;

  std::string sImaDirectory = "";
  unsigned int uTracker = 0;

  // Features
  unsigned int maxTrackedFeatures = 1500;

  // Command options
  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('t', uTracker, "tracker") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--imadir path] \n"
    << "[-t|--tracker Used tracking interface] \n"
    << "\t 0 (default) description based Tracking -> Fast detector + Dipole descriptor\n"
/*
#if defined HAVE_OPENCV
    << "\t 1 image based Tracking -> use OpenCV Pyramidal KLT Tracking\n"
#endif
*/
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  // Call summary
  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImaDirectory << std::endl;

  // Check input folder
  if (sImaDirectory.empty() || !stlplus::is_folder(sImaDirectory))
  {
    std::cerr << "\nIt is an invalid input directory" << std::endl;
    return EXIT_FAILURE;
  }

  // -- Graphics
  if ( !glfwInit() )
  {
    return EXIT_FAILURE;
  }

  CGlWindow window;
  GLuint text2D;

  //---------------------------------------
  // VSSLAM: Visual (SfM_SLAM)
  //---------------------------------------

  // Loop over the image sequence
  //  . detect keypoints
  //  . show them in an openGL window
  //  . track features
  //  . perform VO (and more)

  // Image management
  image::Image<unsigned char> currentImage;
  // Load images from folder
  std::vector<std::string> vec_image = stlplus::folder_files(sImaDirectory);
  // clean invalid image file
  {
    std::vector<std::string> vec_image_;
    for (size_t i = 0; i < vec_image.size(); ++i)
    {
      if (openMVG::image::GetFormat(vec_image[i].c_str()) != openMVG::image::Unknown)
        vec_image_.push_back(vec_image[i]);
    }
    vec_image_.swap(vec_image);
  }
  std::sort(vec_image.begin(), vec_image.end());

  // TODO: Load masks - per camera/image

  // VSSLAM
  using namespace openMVG::VSSLAM;

  // Tracker and Feature detector/matcher interface
  std::unique_ptr<Abstract_Tracker> tracker_ptr;
  std::unique_ptr<Abstract_FeatureExtractor> feat_extractor_ptr;




  return 0;
}
