
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/image/image.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/cameras/cameras.hpp"
#include <software/VO/CGlWindow.hpp>

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include <stdlib.h>
#include <iostream>
#include "openMVG/stl/split.hpp"

// Tracker
#include <openMVG/vsslam/tracking/Tracker_Features.hpp>

#include <openMVG/vsslam/detection/Feat_Extractor_FastDipole.hpp>
#include <openMVG/vsslam/detection/Feat_Extractor_SIFT.hpp>

#include <openMVG/vsslam/matching/Feat_Matcher_Basic.hpp>
#include <openMVG/vsslam/matching/Feat_Matcher_CascadeHashing.hpp>

#include <openMVG/vsslam/SLAM_Monocular.hpp>
#include <openMVG/vsslam/Camera.hpp>
#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include "openMVG/types.hpp"


using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::image;


int n_dummy_param = 0;

/// Check that Kmatrix is a string like "f;0;ppx;0;f;ppy;0;0;1"
/// With f,ppx,ppy as valid numerical value
bool checkIntrinsicStringValidity(const std::string & Kmatrix, double & focal, double & ppx, double & ppy)
{
  std::vector<std::string> vec_str;
  stl::split(Kmatrix, ';', vec_str);
  if (vec_str.size() != 9)  {
    std::cerr << "\n Missing ';' character" << std::endl;
    return false;
  }
  // Check that all K matrix value are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i) {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      std::cerr << "\n Used an invalid not a number character" << std::endl;
      return false;
    }
    if (i==0) focal = readvalue;
    if (i==2) ppx = readvalue;
    if (i==5) ppy = readvalue;
  }
  return true;
}



int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "VISUAL SSLAM -- Tracking demo --" << std::endl;

  CmdLine cmd;

  std::string sImaDirectory = "";
  unsigned int uTracker = 0;

  // Features
  unsigned int maxTrackedFeatures = 1500;

  // Camera data
  std::string sKmatrix;
  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;


  // Command options
  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('t', uTracker, "tracker") );
  cmd.add( make_option('k', sKmatrix, "intrinsics") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--imadir path] \n"
    << "[-t|--tracker Used tracking interface] \n"
    << "\t 0 (default) description based Tracking -> Fast detector + Dipole descriptor\n"
    << "[-k|--intrinsics] Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\"\n"
    << "[-c|--camera_model] Camera model type:\n"
    << "\t 1: Pinhole\n"
    << "\t 2: Pinhole radial 1\n"
    << "\t 3: Pinhole radial 3 (default)\n"
    << "\t 4: Pinhole brown 2\n"
    << "\t 5: Pinhole with a simple Fish-eye distortion\n"
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
            << "--imageDirectory " << sImaDirectory << std::endl
            << "--intrinsics " << sKmatrix << std::endl
            << "--camera_model " << i_User_camera_model << std::endl;

  // Check input folder
  if (sImaDirectory.empty() || !stlplus::is_folder(sImaDirectory))
  {
    std::cerr << "\nIt is an invalid input directory" << std::endl;
    return EXIT_FAILURE;
  }

  // Camera parameters
  openMVG::VSSLAM::CameraParams params_cam_0;
  // cam type
  params_cam_0.camera_model = EINTRINSIC(i_User_camera_model);
  // Basic parameters
  if (sKmatrix.empty() ||
    !checkIntrinsicStringValidity(sKmatrix, params_cam_0.focal, params_cam_0.ppx, params_cam_0.ppy) )
  {
    std::cerr << "\nInvalid K matrix input" << std::endl;
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

  std::cout << "VSSLAM START" << std::endl;
  // Tracker and Feature detector/matcher interface
  std::unique_ptr<Abstract_Tracker> tracker_ptr;
  std::unique_ptr<Abstract_FeatureExtractor> feat_extractor_ptr;
  std::unique_ptr<Abstract_FeatureMatcher> feat_matcher_ptr;
  //feat_matcher_ptr.reset(new Feat_Matcher_Basic());

  switch (uTracker)
  {
    case 0:
      // Set Fast Dipole feature detector/descriptor
      //feat_extractor_ptr.reset(new Feat_Extractor_FastDipole());
      feat_extractor_ptr.reset(new Feat_Extractor_SIFT());
      feat_matcher_ptr.reset(new Feat_Matcher_Basic());
      //feat_matcher_ptr.reset(new Feat_Matcher_CascadeHashing(feat_extractor_ptr.get()));

      tracker_ptr.reset(new Tracker_Features(feat_extractor_ptr.get(),feat_matcher_ptr.get(),maxTrackedFeatures));
      break;
    case 1:
      feat_extractor_ptr.reset(new Feat_Extractor_SIFT());
      feat_matcher_ptr.reset(new Feat_Matcher_CascadeHashing(feat_extractor_ptr.get()));
      //feat_matcher_ptr.reset(new Feat_Matcher_Basic());
      tracker_ptr.reset(new Tracker_Features(feat_extractor_ptr.get(),feat_matcher_ptr.get(),maxTrackedFeatures));
      break;
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
  SLAM_Monocular monocular_slam(tracker_ptr.get());
  monocular_slam.setMapFeatureExtractor(feat_extractor_ptr.get());

  // Load image settings
  if (!params_cam_0.readImageSettings(stlplus::create_filespec( sImaDirectory, vec_image[0] )))
  {
    std::cerr << "\nError reading image header file" << std::endl;
    return EXIT_FAILURE;
  }

  // Create camera model
  if (params_cam_0.checkValidParams())
  {
    params_cam_0.bCalibrated = true;

    if (monocular_slam.createCamera(params_cam_0) < 0)
    {
      std::cerr << "Error: unknown camera model: " << (int) params_cam_0.camera_model << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    std::cerr << "Error: invalid camera parameters"<< std::endl;
    return EXIT_FAILURE;
  }


  if(!monocular_slam.isReady())
  {
    std::cerr << "Error: Monocular SLAM not correctly initialized" << std::endl;
    return EXIT_FAILURE;
  }

  // FOR DISPLAY
  // -----------------
  // -- FRAME BY FRAME processing
  // -----------------
  IndexT frameId = 0;
  IndexT camId = 0;
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

      double startTime = omp_get_wtime();
      // do stuff
      monocular_slam.nextFrame(currentImage, frameId, camId);

      std::cout<<"MonocularSLAM time:"<<omp_get_wtime() - startTime<<"\n";

      Tracker_Features * tr = dynamic_cast<Tracker_Features*>(monocular_slam.tracker_);

      Frame * current_frame = tr->mPrevFrame.get();

      //--
      // Draw feature trajectories
      //--

      IndexT frame_id = tr->mPrevFrame->getFrameId();
      IndexT min_frame_id = frame_id<5 ? 0 : frame_id-5;

      glColor3f(0.f, 1.f, 0.f);
      glLineWidth(2.f);

      if( monocular_slam.tracker_->getTrackingStatus() == Abstract_Tracker::TRACKING_STATUS::OK)
      {
        //std::cout<<"AA\n";
        size_t nn = 0;
        Frame * current_frame = tr->mPrevFrame.get();


        //std::cout<<"AA2\n";
        for (MapLandmark * map_point : current_frame->map_points_)
        {
          if (!map_point)
            continue;

          //std::cout<<"AA3\n";
          LandmarkObservations & map_obs = map_point->obs_;

          //std::cout<<"AA31\n";

          glBegin(GL_LINE_STRIP);


          // Show points with color as degree of connection

          //std::cout<<"AA4\n";
          if (!map_point->isActive()){
            glColor3f(1.f, 0.f, 0.f);
          }
          else{
            switch(map_point->n_all_obs_)
            {
              case 2:
                glColor3f(1.f, 0.f, 0.f);
                break;
              case 3:
                glColor3f(1.f, 1.f, 0.f);
                break;
              case 4:
                glColor3f(0.f, 1.f, 1.f);
                break;
              default:
                glColor3f(0.f, 1.f, 0.f);
                break;
            }
          }
/*
          // How points are associated with current frame
          if (!map_point->isActive())
            continue;

          switch(map_point->association_type_)
          {
            case 1:
              // Initialization point - red
              glColor3f(1.f, 0.f, 0.f);
              break;
            case 2:
              // Motion tracking / reference kf - green
              glColor3f(0.f, 1.f, 0.f);
              break;
            case 3:
              // Map tracking - light blue
              glColor3f(0.f, 1.f, 1.f);
              break;
            case 4:
              // New triangulation - blue
              glColor3f(0.f, 0.f, 1.f);
              break;
            default:
              glColor3f(0.f, 0.f, 0.f);
              break;
          }


*/

          //std::cout<<"AA3\n";

          for(auto obs : map_obs)
          {
            const Vec2 & p0 = (obs.second.frame_ptr->getFeaturePositionDetected(obs.second.feat_id));
            glVertex2f(p0.x(), p0.y());
          }
          glEnd();
          nn++;


        }
        std::cout<<"BB: "<<nn<<"\n";

      }
      for (features::PointFeature pts : current_frame->getRegions()->GetRegionsPositions())
      {
        glPointSize(4.0f);
        glBegin(GL_POINTS);
        glColor3f(1.f, 1.f, 1.f); // Yellow

          glVertex2f(pts.x(), pts.y());
        glEnd();

      }
/*
      for (auto p : monocular_slam.tracker_->mPrevFrame->regions->GetRegionsPositions())
      {
        // draw the current tracked point
        {
          glPointSize(4.0f);
          glBegin(GL_POINTS);
          glColor3f(1.f, 1.f, 0.f); // Yellow
          //const Vec2f & p0 = iter->pos_;
          glVertex2f(p.x(), p.y());
          glEnd();
          //std::cout<<"I: "<<p.x()<<", "<<p.y()<<"\n";
        }
      }*/

      glFlush();
      window.Swap(); // Swap openGL buffer

      std::cout<<"DD\n";
      //if (frameId>0)
      //sleep(2);
    }
  }

  glfwTerminate();

  return 0;
}
