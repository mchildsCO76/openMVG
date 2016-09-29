// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "openMVG/exif/exif_IO_EasyExif.hpp"

#include "openMVG/exif/sensor_width_database/ParseDatabase.hpp"

#include "openMVG/image/image.hpp"
#include "openMVG/stl/split.hpp"

#include "openMVG/sfm/sfm.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <tuple>

#ifdef OPENMVG_USE_CXX11
#include <regex>
#endif

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::image;
using namespace openMVG::sfm;


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


#ifdef OPENMVG_USE_CXX11
typedef tuple <string, EINTRINSIC, std::vector<double> > Regex_Camera_Parameters;
/// Check that sCamsParamsRegex is a string like "imgRegex;camType;Kmatrix"
bool checkCamsRegexStringValidity(const std::string & camsRegex, std::vector<Regex_Camera_Parameters> & vecCamParams)
{
  std::vector<std::string> vec_str;
  stl::split(camsRegex, ';', vec_str);
  
  // Start looping through elements
  for (size_t i = 0; i< vec_str.size();)
  {
    // Go to camera type element (second in the element)
    Regex_Camera_Parameters regex_camParams;
    EINTRINSIC camType;
    std::stringstream ss;
    
    // Check if there is enough parameters for at least one basic camera
    if((vec_str.size()-i)<5)
    {
      std::cerr << "\n Not enough elements in the regex string" << std::endl;
      return false;	  
    }
    // Set regex string to the camera unit
    std::get<0>(regex_camParams) = vec_str[i];
    
    // Read camera type
    i=i+1;    
    // Check if string is integer
    char * p ;
    std::strtol(vec_str[i].c_str(), &p, 10);
    if(*p==0)
    {
      // Convert to string and EINTRINSIC
      int iCamType;
      ss.str(vec_str[i]);
      ss >> iCamType;
      camType = EINTRINSIC(iCamType);
    }
    else{
      std::cerr << "\n Used an invalid not a number character in camera type" << std::endl;
      return false;
    }
    std::get<1>(regex_camParams) = camType;
    
    // Read f,ppx,ppy
    i = i+1;
    double param;
    for(int p_i = 0;p_i<3;p_i++)
    {
      std::get<2>(regex_camParams).push_back(std::stod(vec_str[i+p_i]));
    }

    // Read distortion params
    i = i+3;
    int n_dist_param = 0;
    // Check if there is enough elements in the regex string for distortion parameters
    switch(camType)
    {
      case PINHOLE_CAMERA:
        n_dist_param = 0;
      break;
      case PINHOLE_CAMERA_RADIAL1:
        n_dist_param = 1;
      break;
      case PINHOLE_CAMERA_RADIAL3:
        n_dist_param = 3;
      break;
      case PINHOLE_CAMERA_BROWN:
        n_dist_param = 5;
      break;
      case PINHOLE_CAMERA_FISHEYE:
        n_dist_param = 4;
      break;
      default:
      std::cerr << "Error: unknown camera model: " << (int) camType << std::endl;
      return EXIT_FAILURE;
    }

    // Check if there is sufficient number of parameters
    if((vec_str.size()-i)<n_dist_param)
    {
      std::cerr << "\n Not enough elements in the regex string" << std::endl;
      return false;	  
    }
    // Add all distortion parameters  
    for(int p_i = 0;p_i<n_dist_param;p_i++)
    {
      std::get<2>(regex_camParams).push_back(std::stod(vec_str[i+p_i]));
    }
    i=i+n_dist_param;
    
    // Add the regex unit to the list
    vecCamParams.push_back(regex_camParams);
    
    // Print the summary of the regex camera unit
    std::cout << "\nCamera model: " << (int) camType << std::endl;
    std::cout << "F: " << std::get<2>(regex_camParams).at(0)<<" PPX: "<< std::get<2>(regex_camParams).at(1)<<" PPY: "<< std::get<2>(regex_camParams).at(2)<< std::endl;
    std::cout << "Distortion: ";
    for(int p_i=0;p_i<n_dist_param;p_i++)
    {
      std::cout<<std::get<2>(regex_camParams).at(3+p_i) << " ";
    }
    std::cout << "\nRegex: "<<std::get<0>(regex_camParams) << std::endl;
  }
  return true;
}
#endif

//
// Create the description of an input image dataset for OpenMVG toolsuite
// - Export a SfM_Data file with View & Intrinsic data
//
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImageDir,
    sfileDatabase = "",
    sOutputDir = "",
    sKmatrix;
#ifdef OPENMVG_USE_CXX11
  std::string sCamsParamsRegex;
#endif

  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;

  bool b_Group_camera_model = true;

  double focal_pixels = -1.0;

  cmd.add( make_option('i', sImageDir, "imageDirectory") );
  cmd.add( make_option('d', sfileDatabase, "sensorWidthDatabase") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );
  cmd.add( make_option('f', focal_pixels, "focal") );
  cmd.add( make_option('k', sKmatrix, "intrinsics") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('g', b_Group_camera_model, "group_camera_model") );
#ifdef OPENMVG_USE_CXX11
  cmd.add( make_option('r', sCamsParamsRegex, "regex") );
#endif
  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--imageDirectory]\n"
      << "[-d|--sensorWidthDatabase]\n"
      << "[-o|--outputDirectory]\n"
      << "[-f|--focal] (pixels)\n"
      << "[-k|--intrinsics] Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\"\n"
      << "[-c|--camera_model] Camera model type:\n"
      << "\t 1: Pinhole\n"
      << "\t 2: Pinhole radial 1\n"
      << "\t 3: Pinhole radial 3 (default)\n"
      << "\t 4: Pinhole brown 2\n"
      << "\t 5: Pinhole with a simple Fish-eye distortion\n"
      << "[-g|--group_camera_model]\n"
      << "\t 0-> each view have it's own camera intrinsic parameters,\n"
      << "\t 1-> (default) view can share some camera intrinsic parameters\n"
#ifdef OPENMVG_USE_CXX11
      << "[-r|--regex] Regex: \"{regex;camera_model;f;ppx;ppy;{dist param 1;dist param 2;etc.};]{1,n}\"\n"
#endif
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imageDirectory " << sImageDir << std::endl
            << "--sensorWidthDatabase " << sfileDatabase << std::endl
            << "--outputDirectory " << sOutputDir << std::endl
            << "--focal " << focal_pixels << std::endl
            << "--intrinsics " << sKmatrix << std::endl
            << "--camera_model " << i_User_camera_model << std::endl
            << "--group_camera_model " << b_Group_camera_model << std::endl
#ifdef OPENMVG_USE_CXX11
            << "--regex " << sCamsParamsRegex << std::endl;
#endif

  // Expected properties for each image
  double width = -1, height = -1, focal = -1, ppx = -1,  ppy = -1;

  EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

  if ( !stlplus::folder_exists( sImageDir ) )
  {
    std::cerr << "\nThe input directory doesn't exist" << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutputDir.empty())
  {
    std::cerr << "\nInvalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutputDir ) )
  {
    if ( !stlplus::folder_create( sOutputDir ))
    {
      std::cerr << "\nCannot create output directory" << std::endl;
      return EXIT_FAILURE;
    }
  }
#ifdef OPENMVG_USE_CXX11
  std::vector<Regex_Camera_Parameters> vecCamParams;
  if(sCamsParamsRegex.size()>0){
    if(!checkCamsRegexStringValidity(sCamsParamsRegex,vecCamParams)){
      std::cerr << "\nInvalid Cameras parameters regex string input" << std::endl;
      return EXIT_FAILURE;
    }
  }
#endif

  if (sKmatrix.size() > 0 &&
    !checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy) )
  {
    std::cerr << "\nInvalid K matrix input" << std::endl;
    return EXIT_FAILURE;
  }

  if (sKmatrix.size() > 0 && focal_pixels != -1.0)
  {
    std::cerr << "\nCannot combine -f and -k options" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<Datasheet> vec_database;
  if (!sfileDatabase.empty())
  {
    if ( !parseDatabase( sfileDatabase, vec_database ) )
    {
      std::cerr
       << "\nInvalid input database: " << sfileDatabase
       << ", please specify a valid file." << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::vector<std::string> vec_image = stlplus::folder_files( sImageDir );
  std::sort(vec_image.begin(), vec_image.end());

  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  sfm_data.s_root_path = sImageDir; // Setup main image root_path
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  C_Progress_display my_progress_bar( vec_image.size(),
      std::cout, "\n- Image listing -\n" );
  std::ostringstream error_report_stream;
  for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin();
    iter_image != vec_image.end();
    ++iter_image, ++my_progress_bar )
  {
    // Read meta data to fill camera parameter (w,h,focal,ppx,ppy) fields.
    width = height = ppx = ppy = focal = -1.0;

    const std::string sImageFilename = stlplus::create_filespec( sImageDir, *iter_image );
    const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

    // Test if the image format is supported:
    if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown)
    {
      error_report_stream
          << sImFilenamePart << ": Unkown image file format." << "\n";
      continue; // image cannot be opened
    }

    if(sImFilenamePart.find("mask.png") != std::string::npos
       || sImFilenamePart.find("_mask.png") != std::string::npos)
    {
      error_report_stream
          << sImFilenamePart << " is a mask image" << "\n";
      continue;
    }

    ImageHeader imgHeader;
    if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
      continue; // image cannot be read

    width = imgHeader.width;
    height = imgHeader.height;
    ppx = width / 2.0;
    ppy = height / 2.0;

    std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
    exifReader->open( sImageFilename );

    const bool bHaveValidExifMetadata =
      exifReader->doesHaveExifInfo()
      && !exifReader->getModel().empty();

#ifdef OPENMVG_USE_CXX11
    bool bImageMatched = false;
    // Check if image fits any regex patterns
    std::vector<double> cam_image_param;
    if(vecCamParams.size()>0)
    {
      // Loop through regex patterns
      for(std::vector<Regex_Camera_Parameters>::const_iterator iter_cParam = vecCamParams.begin();iter_cParam!=vecCamParams.end(); ++iter_cParam)
      {
        const std::regex imgNameRegex(std::get<0>(*iter_cParam),std::regex_constants::extended);
        std::smatch smBaseCamMatch;
        // Check if it match the regex
        if (std::regex_match(*iter_image, smBaseCamMatch, imgNameRegex))
        {
          e_User_camera_model = EINTRINSIC(std::get<1>(*iter_cParam));
          // Get intrinsic parameters of the model
          cam_image_param = std::get<2>(*iter_cParam);
          bImageMatched = true;
          break;
        }
      }
    }
#endif

    // Build intrinsic parameter related to the view
    std::shared_ptr<IntrinsicBase> intrinsic (NULL);

#ifdef OPENMVG_USE_CXX11    
    // Image has not been matched with regex
    if(!bImageMatched)
    {
#endif
      // Consider the case where the focal is provided manually
      if ( !bHaveValidExifMetadata || focal_pixels != -1)
      {
        if (sKmatrix.size() > 0) // Known user calibration K matrix
        {
          if (!checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy))
            focal = -1.0;
        }
        else // User provided focal length value
          if (focal_pixels != -1 )
            focal = focal_pixels;
      }
      else // If image contains meta data
      {
        const std::string sCamModel = exifReader->getModel();

        // Handle case where focal length is equal to 0
        if (exifReader->getFocal() == 0.0f)
        {
          error_report_stream
            << stlplus::basename_part(sImageFilename) << ": Focal length is missing." << "\n";
          focal = -1.0;
        }
        else
        // Create the image entry in the list file
        {
          Datasheet datasheet;
          if ( getInfo( sCamModel, vec_database, datasheet ))
          {
            // The camera model was found in the database so we can compute it's approximated focal length
            const double ccdw = datasheet.sensorSize_;
            focal = std::max ( width, height ) * exifReader->getFocal() / ccdw;
          }
          else
          {
            error_report_stream
              << stlplus::basename_part(sImageFilename)
              << "\" model \"" << sCamModel << "\" doesn't exist in the database" << "\n"
              << "Please consider add your camera model and sensor width in the database." << "\n";
          }
        }
      }
      e_User_camera_model = EINTRINSIC(i_User_camera_model);
      
      
      if (focal > 0 && ppx > 0 && ppy > 0 && width > 0 && height > 0)
      {
        // Create the desired camera type
        switch(e_User_camera_model)
        {
          case PINHOLE_CAMERA:
          intrinsic = std::make_shared<Pinhole_Intrinsic>
              (width, height, focal, ppx, ppy);
          break;
          case PINHOLE_CAMERA_RADIAL1:
            intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
              (width, height, focal, ppx, ppy, 0.0); // setup no distortion as initial guess
          break;
          case PINHOLE_CAMERA_RADIAL3:
            intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
              (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0);  // setup no distortion as initial guess
          break;
          case PINHOLE_CAMERA_BROWN:
            intrinsic =std::make_shared<Pinhole_Intrinsic_Brown_T2>
              (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
          break;
          case PINHOLE_CAMERA_FISHEYE:
            intrinsic =std::make_shared<Pinhole_Intrinsic_Fisheye>
              (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
          break;
          default:
            std::cerr << "Error: unknown camera model: " << (int) e_User_camera_model << std::endl;
            return EXIT_FAILURE;
        }
      }
#ifdef OPENMVG_USE_CXX11
    }
    else
    {
      focal = cam_image_param.at(0);
      ppx = cam_image_param.at(1);
      ppy = cam_image_param.at(2);
      if (focal > 0 && ppx > 0 && ppy > 0 && width > 0 && height > 0)
      {
        // Create the desired camera type
        switch(e_User_camera_model)
        {
          case PINHOLE_CAMERA:
          intrinsic = std::make_shared<Pinhole_Intrinsic>
              (width, height, focal, ppx, ppy);
          break;
          case PINHOLE_CAMERA_RADIAL1:
            intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
              (width, height, focal, ppx, ppy, cam_image_param.at(3)); // setup no distortion as initial guess
          break;
          case PINHOLE_CAMERA_RADIAL3:
            intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
              (width, height, focal, ppx, ppy, cam_image_param.at(3), cam_image_param.at(4), cam_image_param.at(5));  // setup no distortion as initial guess
          break;
          case PINHOLE_CAMERA_BROWN:
            intrinsic =std::make_shared<Pinhole_Intrinsic_Brown_T2>
              (width, height, focal, ppx, ppy, cam_image_param.at(3), cam_image_param.at(4), cam_image_param.at(5), cam_image_param.at(6), cam_image_param.at(7)); // setup no distortion as initial guess
          break;
          case PINHOLE_CAMERA_FISHEYE:
            intrinsic =std::make_shared<Pinhole_Intrinsic_Fisheye>
              (width, height, focal, ppx, ppy, cam_image_param.at(3), cam_image_param.at(4), cam_image_param.at(5), cam_image_param.at(6)); // setup no distortion as initial guess
          break;
          default:
            std::cerr << "Error: unknown camera model: " << (int) e_User_camera_model << std::endl;
            return EXIT_FAILURE;
        }
      }
    }
#endif

    // Build the view corresponding to the image
    View v(*iter_image, views.size(), views.size(), views.size(), width, height);

    // Add intrinsic related to the image (if any)
    if (intrinsic == NULL)
    {
      //Since the view have invalid intrinsic data
      // (export the view, with an invalid intrinsic field value)
      v.id_intrinsic = UndefinedIndexT;
    }
    else
    {
      // Add the defined intrinsic to the sfm_container
      intrinsics[v.id_intrinsic] = intrinsic;
    }

    // Add the view to the sfm_container
    views[v.id_view] = std::make_shared<View>(v);
  }

  // Display saved warning & error messages if any.
  if (!error_report_stream.str().empty())
  {
    std::cerr
      << "\nWarning & Error messages:" << std::endl
      << error_report_stream.str() << std::endl;
  }

  // Group camera that share common properties if desired (leads to more faster & stable BA).
  if (b_Group_camera_model)
  {
    GroupSharedIntrinsics(sfm_data);
  }

  // Store SfM_Data views & intrinsic data
  if (!Save(
    sfm_data,
    stlplus::create_filespec( sOutputDir, "sfm_data.json" ).c_str(),
    ESfM_Data(VIEWS|INTRINSICS)))
  {
    return EXIT_FAILURE;
  }

  std::cout << std::endl
    << "SfMInit_ImageListing report:\n"
    << "listed #File(s): " << vec_image.size() << "\n"
    << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << std::endl;

  return EXIT_SUCCESS;
}
