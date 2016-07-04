

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/cameras/cameras.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"


using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "openMVG/numeric/numeric.h"


int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";
  bool bEstimateLandmarks = false;
  std::string sOutDir = "";

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
  cmd.add( make_option('l', bEstimateLandmarks, "estimateLandmarks") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--sfmdata] filename, the SfM_Data file to use\n"
      << "[-f|--refineIntrinsics] Intrinsic parameters refinement option\n"
        << "\t ADJUST_ALL -> refine all existing parameters (default) \n"
        << "\t NONE -> intrinsic parameters are held as constant\n"
        << "\t ADJUST_FOCAL_LENGTH -> refine only the focal length\n"
        << "\t ADJUST_PRINCIPAL_POINT -> refine only the principal point position\n"
        << "\t ADJUST_DISTORTION -> refine only the distortion coefficient(s) (if any)\n"
        << "\t -> NOTE: options can be combined thanks to '|'\n"
        << "\t ADJUST_FOCAL_LENGTH|ADJUST_PRINCIPAL_POINT\n"
        <<      "\t\t-> refine the focal length & the principal point position\n"
        << "\t ADJUST_FOCAL_LENGTH|ADJUST_DISTORTION\n"
        <<      "\t\t-> refine the focal length & the distortion coefficient(s) (if any)\n"
        << "\t ADJUST_PRINCIPAL_POINT|ADJUST_DISTORTION\n"
        <<      "\t\t-> refine the principal point position & the distortion coefficient(s) (if any)\n"
      << "[-l|--estimateLandmarks] estimate uncertainty of landmarks\n"
      << "[-o|--outdir] path\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
    cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--sfmdata " << sSfM_Data_Filename << std::endl
            << "--outdir " << sOutDir << std::endl;

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  // Read the SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "\n\n-------------------------------" << "\n"
      << "Estimate Uncertainty"<<"\n"
      << "-------------------------------" << "\n";
  Bundle_Adjustment_Ceres::BA_Ceres_options options;
  Bundle_Adjustment_Ceres bundle_adjustment_obj(options);
  const Optimize_Options ba_refine_options
    ( intrinsic_refinement_options,
      Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
      Structure_Parameter_Type::ADJUST_ALL // Adjust scene structure
    );
  bundle_adjustment_obj.EstimateUncertainty(sfm_data, ba_refine_options,bEstimateLandmarks);


  //-- Export to disk computed scene (data & visualizable results)
  std::cout << "...Export SfM_Data to disk." << std::endl;
  Save(sfm_data,
    stlplus::create_filespec(sOutDir, "sfm_data", ".json"),
    ESfM_Data(ALL|UNCERTAINTIES));


  return EXIT_SUCCESS;
}
