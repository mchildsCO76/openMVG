
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

#include <openMVG/vsslam/matching/Feat_Matcher_CascadeHashing.hpp>

#include <openMVG/vsslam/SLAM_Monocular.hpp>
#include <openMVG/vsslam/Camera.hpp>
#include <openMVG/vsslam/Frame.hpp>
#include <openMVG/vsslam/VSSLAM_Data.hpp>
#include "openMVG/types.hpp"

#include <openMVG/vsslam/Export.hpp>
#include <openMVG/vsslam/matching/Feat_Matcher_Basic.hpp>


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

void DrawCircle(float cx, float cy, float r)
{
    int num_segments = 20;
    glBegin(GL_LINE_LOOP);
    for(int ii = 0; ii < num_segments; ii++)
    {
        float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle

        float x = r * cosf(theta);//calculate the x component
        float y = r * sinf(theta);//calculate the y component

        glVertex2f(x + cx, y + cy);//output vertex

    }
    glEnd();
}


void DrawSteps()
{


}



int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "VISUAL SSLAM -- Tracking demo --" << std::endl;

  CmdLine cmd;

  std::string sImaDirectory = "";
  std::string sImaMask = "";
  unsigned int uTracker = 0;

  // Features
  unsigned int maxTrackedFeatures = 1500;

  // Camera data
  std::string sKmatrix;
  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;


  // Command options
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
            << "--imageMask " << sImaMask << std::endl
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
  bool b_use_mask = false;
  image::Image<unsigned char> maskImage;
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


  // TODO: Load masks - per camera/image
  if (!sImaMask.empty() && stlplus::file_exists(sImaMask))
  {
    if (openMVG::image::GetFormat(sImaMask.c_str()) == openMVG::image::Unknown)
      std::cout << "\nMask image path is invalid! Not using mask image!" << std::endl;
    else
    {
      if (!(openMVG::image::ReadImage( sImaMask.c_str(), &maskImage)))
      {
        std::cout << "\nMask image is invalid! Not using mask image!" << std::endl;
      }
      else
      {
        std::cout << "\nUsing mask image: "<<sImaMask<<"!\n" << std::endl;
        b_use_mask = true;
      }
    }

  }



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
      //feat_matcher_ptr.reset(new Feat_Matcher_CascadeHashing(feat_extractor_ptr.get()));
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

    if (monocular_slam.createCamera(params_cam_0, b_use_mask ? &maskImage : nullptr) < 0)
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
  std::shared_ptr<Frame> display_prev_frame;
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
      Abstract_Tracker::TRACKING_STATUS pre_status = monocular_slam.tracker_->getTrackingStatus();

      double startTime = omp_get_wtime();
      // do stuff
      monocular_slam.nextFrame(currentImage, frameId, camId);

      std::cout<<"MonocularSLAM time:"<<omp_get_wtime() - startTime<<"\n";

      Tracker_Features * tr = dynamic_cast<Tracker_Features*>(monocular_slam.tracker_);

      Frame * display_current_frame = tr->mPrevFrame.get();

      //--
      // Draw feature trajectories
      //--

      IndexT frame_id = tr->mPrevFrame->getFrameId();
      IndexT min_frame_id = frame_id<5 ? 0 : frame_id-5;

      glColor3f(0.f, 1.f, 0.f);
      glLineWidth(2.f);


      // Draw all features in current frame - YELLOW
      for (size_t map_point_i = 0; map_point_i<display_current_frame->map_points_.size(); map_point_i++)
      {
        glPointSize(2.0f);
        glColor3f(1.f, 1.f, 0.f);
        glBegin(GL_POINTS);
        Vec2 pt = display_current_frame->getFeaturePosition(map_point_i);
        glVertex2f(pt(0),pt(1));
        glEnd();
      }


      if (false)
      {
        if( pre_status == Abstract_Tracker::TRACKING_STATUS::NOT_INIT || pre_status == Abstract_Tracker::TRACKING_STATUS::INIT)
        {
          ////////////////////////////////////////////////////////////////////////
          // Show initial matches &&
          ////////////////////////////////////////////////////////////////////////
          size_t step_i = 0;
          for (step_i = 0; step_i < monocular_slam.tracker_->display_pt2d_A.size(); step_i++)
          {
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

            // Draw all features in current frame - YELLOW
            for (size_t map_point_i = 0; map_point_i<display_current_frame->map_points_.size(); map_point_i++)
            {
              glPointSize(2.0f);
              glColor3f(1.f, 1.f, 0.f);
              glBegin(GL_POINTS);
              Vec2 pt = display_current_frame->getFeaturePosition(map_point_i);
              glVertex2f(pt(0),pt(1));
              glEnd();
            }

            for (size_t i=0; i<monocular_slam.tracker_->display_pt2d_A[step_i].size(); ++i)
            {
              Vec2 ptA = monocular_slam.tracker_->display_pt2d_A[step_i][i];
              Vec2 ptB = monocular_slam.tracker_->display_pt2d_B[step_i][i];
              Vec2 ptC = monocular_slam.tracker_->display_pt2d_C[step_i][i];


              glPointSize(4.0f);
              glBegin(GL_POINTS);

              // feature matched in prev frame  - RED
              glColor3f(1.f, 0.f, 0.f);
              glVertex2f(ptA(0),ptA(1));

              glEnd();

              DrawCircle(ptA(0),ptA(1),ptC(0));

              glBegin(GL_POINTS);
              // Projected point to current frame  - BLUE
              glColor3f(0.f, 1.f, 0.f);
              glVertex2f(ptB(0),ptB(1));

              glEnd();

              DrawCircle(ptB(0),ptB(1),ptC(1));

              // Line from point in prev to predicted - WHITE
              glColor3f(1.f, 1.f, 1.f);
              glBegin(GL_LINE_STRIP);
              glVertex2f(ptA(0),ptA(1));
              glVertex2f(ptB(0),ptB(1));
              glEnd();

            }

            glFlush();
            window.Swap(); // Swap openGL buffer
            if (step_i<(monocular_slam.tracker_->display_pt2d_A.size()-1))
              sleep(3);
          }
        }
        else if ( pre_status == Abstract_Tracker::TRACKING_STATUS::OK)
        {

          ////////////////////////////////////////////////////////////////////////
          // Show putative matches of matching with MM
          ////////////////////////////////////////////////////////////////////////
          size_t step_i;
          size_t first_part_i =  monocular_slam.tracker_->display_iterations[0] + monocular_slam.tracker_->display_iterations[1] + monocular_slam.tracker_->display_iterations[2] + monocular_slam.tracker_->display_iterations[3];
          for (step_i = 0; step_i <first_part_i; step_i++)
          {
            std::cout<<"Step :"<<step_i<<"----------------\n";
            std::cout<<" Title: "<<monocular_slam.tracker_->display_text[step_i]<<"\n";
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

            // Draw all features in current frame - YELLOW
            for (size_t map_point_i = 0; map_point_i<display_current_frame->map_points_.size(); map_point_i++)
            {
              glPointSize(2.0f);
              glColor3f(1.f, 1.f, 0.f);
              glBegin(GL_POINTS);
              Vec2 pt = display_current_frame->getFeaturePosition(map_point_i);
              glVertex2f(pt(0),pt(1));
              glEnd();
            }

            for (size_t i=0; i<monocular_slam.tracker_->display_pt2d_A[step_i].size(); ++i)
            {
              Vec2 ptA = monocular_slam.tracker_->display_pt2d_A[step_i][i];
              Vec2 ptB = monocular_slam.tracker_->display_pt2d_B[step_i][i];
              Vec2 ptC = monocular_slam.tracker_->display_pt2d_C[step_i][i];
              float pt_s = monocular_slam.tracker_->display_size_A[step_i][i];


              glPointSize(4.0f);
              glBegin(GL_POINTS);

              // Point in prev frame  - RED
              glColor3f(1.f, 0.f, 0.f);
              if (ptB(0)!=-1 && ptB(1)!=-1)
              {
                glVertex2f(ptB(0),ptB(1));
              }
              // Matched point in current frame - BLUE
              if (ptC(0)!=-1 && ptC(1)!=-1)
              {
                glColor3f(0.f, 0.f, 1.f);
                glVertex2f(ptC(0),ptC(1));
              }

              // Predicted location of the point based on prev_frame and camera pose  - GREEN
              if (ptC(0)!=-1 && ptC(1)!=-1)
              {
                glColor3f(0.f, 1.f, 0.f);
              }
              else
              {
                glColor3f(1.f, 0.f, 0.f);
              }
              glVertex2f(ptA(0),ptA(1));

              glEnd();

              // Line from predicted to matched - green
              if (ptC(0)!=-1 && ptC(1)!=-1)
              {
                glColor3f(0.f, 1.f, 0.f);
                glBegin(GL_LINE_STRIP);
                glVertex2f(ptA(0),ptA(1));
                glVertex2f(ptC(0),ptC(1));
                glEnd();
              }
              // Area of search around predicted point
              if (ptC(0)!=-1 && ptC(1)!=-1)
              {
                glColor3f(0.f, 1.f, 0.f);
              }
              else
              {
                glColor3f(1.f, 0.f, 0.f);
              }
              DrawCircle(ptA(0),ptA(1),pt_s);
            }

            glFlush();
            window.Swap(); // Swap openGL buffer

              sleep(3);
          }

          ////////////////////////////////////////////////////////////////////////
          // Show matches with epipolar lines
          ////////////////////////////////////////////////////////////////////////
          size_t second_part_i = first_part_i + monocular_slam.tracker_->display_iterations[4];
          for (step_i = first_part_i; step_i < second_part_i; step_i++)
          {
            std::cout<<"Step :"<<step_i<<"----------------\n";
            std::cout<<" Title: "<<monocular_slam.tracker_->display_text[step_i]<<"\n";
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

            // Draw all features in current frame - YELLOW
            for (size_t map_point_i = 0; map_point_i<display_current_frame->map_points_.size(); map_point_i++)
            {
              glPointSize(2.0f);
              glColor3f(1.f, 1.f, 0.f);
              glBegin(GL_POINTS);
              Vec2 pt = display_current_frame->getFeaturePosition(map_point_i);
              glVertex2f(pt(0),pt(1));
              glEnd();
            }

            glPointSize(8.0f);
            glColor3f(0.f, 1.f, 1.f);
            glBegin(GL_POINTS);
            Vec2 pt = monocular_slam.tracker_->display_pt2d_A[step_i][0];
            if (pt(0)>0 && pt(0)<currentImage.Width() && pt(1)>0 && pt(1)<currentImage.Height())
              glVertex2f(pt(0),pt(1));
            glEnd();

            for (size_t i=1; i<monocular_slam.tracker_->display_pt2d_A[step_i].size(); ++i)
            {
              Vec2 ptA = monocular_slam.tracker_->display_pt2d_A[step_i][i];
              Vec2 ptB = monocular_slam.tracker_->display_pt2d_B[step_i][i];
              float pt_s = monocular_slam.tracker_->display_size_A[step_i][i];


              glPointSize(4.0f);
              glBegin(GL_POINTS);

              // Point in prev frame  - RED
              glColor3f(1.f, 0.f, 0.f);
              glVertex2f(ptB(0),ptB(1));

              // Matched point in current frame - BLUE
              glColor3f(0.f, 0.f, 1.f);
              glVertex2f(ptA(0),ptA(1));

              glEnd();

              // Line from point in prev to predicted - WHITE
              glColor3f(1.f, 1.f, 1.f);
              glBegin(GL_LINE_STRIP);
              glVertex2f(ptA(0),ptA(1));
              glVertex2f(ptB(0),ptB(1));
              glEnd();

            }

            glFlush();
            window.Swap(); // Swap openGL buffer

            sleep(3);
          }




          ////////////////////////////////////////////////////////////////////////
          // Show putative new triangulated pts
          ////////////////////////////////////////////////////////////////////////

          size_t third_part_i = second_part_i + monocular_slam.tracker_->display_iterations[5];
          for (step_i = second_part_i; step_i < third_part_i; step_i++)
          {
            std::cout<<"Step :"<<step_i<<"----------------\n";
            std::cout<<" Title: "<<monocular_slam.tracker_->display_text[step_i]<<"\n";
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

            // Draw all features in current frame - YELLOW
            for (size_t map_point_i = 0; map_point_i<display_current_frame->map_points_.size(); map_point_i++)
            {
              glPointSize(2.0f);
              glColor3f(1.f, 1.f, 0.f);
              glBegin(GL_POINTS);
              Vec2 pt = display_current_frame->getFeaturePosition(map_point_i);
              glVertex2f(pt(0),pt(1));
              glEnd();
            }

            for (size_t i=0; i<monocular_slam.tracker_->display_pt2d_A[step_i].size(); ++i)
            {
              Vec2 ptA = monocular_slam.tracker_->display_pt2d_A[step_i][i];
              Vec2 ptC = monocular_slam.tracker_->display_pt2d_C[step_i][i];
              float pt_s = monocular_slam.tracker_->display_size_A[step_i][i];


              glPointSize(4.0f);
              glBegin(GL_POINTS);


              // New triangulated pts
              if (ptC(0) == 6)
                glColor3f(0.f, 0.f, 1.f); // 2 obs BLUE
              else if (ptC(0) == 7)
                glColor3f(0.f, 1.f, 0.f); // more obs GREEN
              else
                glColor3f(1.f, 0.f, 0.f);  // error RED

              glVertex2f(ptA(0),ptA(1));

              glEnd();

              // Area of search around predicted point
              DrawCircle(ptA(0),ptA(1),pt_s);
            }

            glFlush();
            window.Swap(); // Swap openGL buffer

            if (step_i < monocular_slam.tracker_->display_pt2d_A.size()-1)
              sleep(3);
          }
        }
      }

      monocular_slam.tracker_->display_pt2d_A.clear();
      monocular_slam.tracker_->display_pt2d_B.clear();
      monocular_slam.tracker_->display_pt2d_C.clear();
      monocular_slam.tracker_->display_size_A.clear();
      monocular_slam.tracker_->display_text.clear();

      if (true)
      {
      bool b_snake_trail = true;

      if( monocular_slam.tracker_->getTrackingStatus() == Abstract_Tracker::TRACKING_STATUS::OK)
      {
        for (size_t map_point_i = 0; map_point_i<display_current_frame->map_points_.size(); map_point_i++)
        {
          MapLandmark * map_point = display_current_frame->map_points_[map_point_i];

          if (!map_point)
            continue;

          // Define color based on the type of association of the point
          Vec3 color;
          if (!map_point->isActive())
          {
            color = Vec3(1.0,1.0,1.0);  // WHITE
          }
          else
          {
            switch(map_point->association_type_)
            {
              case 1:
              // Initialization points
                color = Vec3(1.0,0.0,0.0);  // RED
              break;
              case 2:
              // Motion Model points
                color = Vec3(0.0,1.0,0.0);  // GREEN
              break;
              case 3:
              // reference frame points
                color = Vec3(1.0,1.0,0.0);  // YELLOW
              break;
              case 4:
              // Reference map points
                color = Vec3(1.0,0.0,1.0);  // MAGENTA
              break;
              case 5:
              // Local map points
                color = Vec3(0.0,1.0,1.0);  // CYAN
              break;
              case 6:
              // Triangulated 2 points
                color = Vec3(1.0,1.0,1.0);  // WHITE
              break;
              case 7:
              // Triangulated multiple points
                color = Vec3(0.0,0.0,1.0);  // BLUE
              break;
            }
          }
          glColor3f(color(0),color(1),color(2));


          // Get projection of the point
          Vec3 pt_3D_frame;
          PoseEstimator::getRelativePointPosition(map_point->X_,map_point->ref_frame_,pt_3D_frame,display_current_frame);

          IntrinsicBase * cam_intrinsic = display_current_frame->getCameraIntrinsics();
          Vec2 pt_3D_frame_projected = cam_intrinsic->cam2ima(pt_3D_frame.hnormalized());

          // Draw point
          glPointSize(4.0f);
          glBegin(GL_POINTS);
          glVertex2f(pt_3D_frame_projected(0),pt_3D_frame_projected(1));
          glEnd();

          // Draw Scale
          //DrawCircle(pt_3D_frame_projected(0),pt_3D_frame_projected(1),display_current_frame->getFeatureScale(map_point_i));
          DrawCircle(pt_3D_frame_projected(0),pt_3D_frame_projected(1),4);

          if (b_snake_trail)
          {
            // Show track trail through all previous tracks
            LandmarkObservations & map_obs = map_point->obs_;

            // Show track trail through all previous tracks
            glBegin(GL_LINE_STRIP);
            for(auto obs : map_obs)
            {
              const Vec2 & p0 = (obs.second.frame_ptr->getFeaturePositionDetected(obs.second.feat_id));
              // Check if projection is actually in the image borders
              if (!display_current_frame->isPointInFrame(p0))
                continue;
              glVertex2f(p0.x(), p0.y());
            }
            glEnd();


            glPointSize(2.0f);
            glColor3f(1.0,1.0,1.0);
            glBegin(GL_POINTS);
            for(auto obs : map_obs)
            {
              const Vec2 & p0 = (obs.second.frame_ptr->getFeaturePositionDetected(obs.second.feat_id));
              // Check if projection is actually in the image borders
              if (!display_current_frame->isPointInFrame(p0))
                continue;
              glVertex2f(p0.x(), p0.y());
            }
            glEnd();
          }
        }
      }
      }
      display_prev_frame = monocular_slam.tracker_->mPrevFrame->share_ptr();
      glFlush();
      window.Swap(); // Swap openGL buffer

      //if (frameId>0)
      //sleep(1);
      std::cout<<"Press ENTER to continue....."<<std::endl<<std::endl;
      std::cin.ignore(1);
    }

  }

  glfwTerminate();

  return 0;
}
