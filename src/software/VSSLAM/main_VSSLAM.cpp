
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
#include <software/VSSLAM/slam/Tracker_Features.hpp>
#include <software/VSSLAM/slam/Feat_Extractor_FastDipole.hpp>

#include <software/VSSLAM/slam/SLAM_Monocular.hpp>


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

  switch (uTracker)
  {
    case 0:
      tracker_ptr.reset(new Tracker_Features);
      // Set Fast Dipole feature detector/descriptor
      feat_extractor_ptr.reset(new Feat_Extractor_FastDipole);
      dynamic_cast<Tracker_Features*>(tracker_ptr.get())->setFeatureExtractor(feat_extractor_ptr.get());
    break;
    default:
    std::cerr << "Unknow tracking method" << std::endl;
    return EXIT_FAILURE;
  }

  if (!tracker_ptr)
  {
    std::cerr << "Cannot instantiate the tracking interface" << std::endl;
    return EXIT_FAILURE;
  }


  // Initialize the monocular tracking framework
  SLAM_Monocular monocular_slam(tracker_ptr.get(), maxTrackedFeatures);

  size_t frameId = 0;
  for (std::vector<std::string>::const_iterator iterFile = vec_image.begin();
    iterFile != vec_image.end(); ++iterFile, ++frameId)
  {
    const std::string sImageFilename = stlplus::create_filespec( sImaDirectory, *iterFile );
    if (openMVG::image::ReadImage( sImageFilename.c_str(), &currentImage))
    {
      if (window._height < 0)
      {
        // no window created yet, initialize it with the first frame

        const double aspect_ratio = currentImage.Width()/(double)currentImage.Height();
        window.Init(640, 640/aspect_ratio, "VisualOdometry--TrackingViewer");
        glGenTextures(1,&text2D);             //allocate the memory for texture
        glBindTexture(GL_TEXTURE_2D, text2D); //Binding the texture
        glEnable(GL_TEXTURE_2D);              //Enable texture
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE,
          currentImage.Width(), currentImage.Height(), 0,
          GL_LUMINANCE, GL_UNSIGNED_BYTE, currentImage.data());
      }
      glBindTexture(GL_TEXTURE_2D, text2D); //Binding the texture
      glEnable(GL_TEXTURE_2D);              //Enable texture
      //-- Update the openGL texture with the current frame pixel values
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0,
        currentImage.Width(), currentImage.Height(),
        GL_LUMINANCE, GL_UNSIGNED_BYTE,
        currentImage.data());

      //-- Draw the current image
      window.SetOrtho(currentImage.Width(), currentImage.Height());
      window.DrawFullScreenTexQuad(currentImage.Width(), currentImage.Height());
      glDisable(GL_TEXTURE_2D);

      // Clear the depth buffer so the drawn image become the background
      glClear(GL_DEPTH_BUFFER_BIT);
      glDisable(GL_LIGHTING);

      //--
      //-- Feature tracking
      //    . track features
      //    . if some tracks are cut, detect and insert new features
      //--

      monocular_slam.nextFrame(currentImage, frameId);

      //--
      // Draw feature trajectories
      //--
      glColor3f(0.f, 1.f, 0.f);
      glLineWidth(2.f);

      if(monocular_slam.tracker_->init_ref_frame && monocular_slam.tracker_->mPrevFrame){
      for (auto track : monocular_slam.tracker_->feat_cur_prev_matches_ids)
      {
        // Draw the line from prev
        glBegin(GL_LINE_STRIP);
        glColor3f(0.f, 1.f, 0.f);
        const Vec2 & p0 = monocular_slam.tracker_->init_ref_frame->regions->GetRegionPosition(track.second);
        const Vec2 & p1 = monocular_slam.tracker_->mPrevFrame->regions->GetRegionPosition(track.first);
        std::cout<<"A: "<<p0(0)<< ", "<<p0(1)<<" :: "<<p1(0)<< ", "<<p1(1)<<"\n";
        glVertex2f(p0(0), p0(1));
        glVertex2f(p1(0), p1(1));
        glEnd();
      }
      }
      else
      {
        std::cout<<"AABAB\n";
      }
/*
      for (auto p : monocular_slam.tracker_->mCurrentFrame->regions->GetRegionsPositions())
      {
        // draw the current tracked point
        {
          glPointSize(4.0f);
          glBegin(GL_POINTS);
          glColor3f(1.f, 1.f, 0.f); // Yellow
          //const Vec2f & p0 = iter->pos_;
          glVertex2f(p.x(), p.y());
          glEnd();
        }
      }
*/
      glFlush();
      window.Swap(); // Swap openGL buffer
      sleep(2);
    }
  }

  glfwTerminate();

  return 0;
}
