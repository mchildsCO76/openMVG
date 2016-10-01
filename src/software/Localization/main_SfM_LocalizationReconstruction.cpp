
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <openMVG/sfm/sfm.hpp>
#include <openMVG/features/features.hpp>
#include <nonFree/sift/SIFT_describer.hpp>
#include <openMVG/image/image.hpp>
#include <software/SfM/SfMPlyHelper.hpp>

#include <openMVG/system/timer.hpp>
#include "openMVG/stl/stl.hpp"
//- Robust estimation - LMeds (since no threshold can be defined)
#include "openMVG/robust_estimation/robust_estimator_LMeds.hpp"
#include "openMVG/geometry/Similarity3_Kernel.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"
#include "openMVG/sfm/sfm_data_filters.hpp"

using namespace openMVG;
using namespace openMVG::sfm;
using namespace openMVG::cameras;
using namespace openMVG::geometry;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>

//  Adjust camera poses (ones that are not used in the localization)  
bool AdjustReconstructionWithSomeFixedCameras
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  const cameras::Intrinsic_Parameter_Type &intrinsic_refinement_options,
  const Hash_Map<IndexT, Pose3> &poses_estimated
)
{ 
  // BA Options
  Bundle_Adjustment_Ceres::BA_Ceres_options ceres_options;
  if ( sfm_data.GetPoses().size() > 100 &&
      (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE) ||
       ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE) ||
       ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
      )
  // Enable sparse BA only if a sparse lib is available and if there more than 100 poses
  {
    ceres_options.preconditioner_type_ = ceres::JACOBI;
    ceres_options.linear_solver_type_ = ceres::SPARSE_SCHUR;
  }
  else
  {
    ceres_options.linear_solver_type_ = ceres::DENSE_SCHUR;
  }
  
  const Optimize_Options options
    ( intrinsic_refinement_options,
      Extrinsic_Parameter_Type::ADJUST_ALL, // Adjust camera motion
      Structure_Parameter_Type::ADJUST_ALL, // Adjust scene structure
      Control_Point_Parameter(),
      false
    );
    
  
  
  
  
  
    //----------
  // Add camera parameters
  // - intrinsics
  // - poses [R|t]

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  //----------
  ceres::Problem problem;
  
  // Data wrapper for refinement:
  Hash_Map<IndexT, std::vector<double> > map_intrinsics;
  Hash_Map<IndexT, std::vector<double> > map_poses;

  // Setup Poses data & subparametrization
  for (const auto & pose_it : sfm_data.poses)
  {
    const IndexT indexPose = pose_it.first;

    const Pose3 & pose = pose_it.second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE || poses_estimated.count(indexPose) != 0)
    {
      // set the whole parameter block as constant for best performance
      problem.SetParameterBlockConstant(parameter_block);
    }
    else  // Subset parametrization
    {
      std::vector<int> vec_constant_extrinsic;
      // If we adjust only the translation, we must set ROTATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_TRANSLATION)
      {
        // Subset rotation parametrization
        vec_constant_extrinsic.push_back(0);
        vec_constant_extrinsic.push_back(1);
        vec_constant_extrinsic.push_back(2);
      }
      // If we adjust only the rotation, we must set TRANSLATION as constant
      if (options.extrinsics_opt == Extrinsic_Parameter_Type::ADJUST_ROTATION)
      {
        // Subset translation parametrization
        vec_constant_extrinsic.push_back(3);
        vec_constant_extrinsic.push_back(4);
        vec_constant_extrinsic.push_back(5);
      }
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(parameter_block, subset_parameterization);
      }
    }
  }

  // Setup Intrinsics data & subparametrization
  for (const auto & intrinsic_it : sfm_data.intrinsics)
  {
    const IndexT indexCam = intrinsic_it.first;

    if (isValid(intrinsic_it.second->getType()))
    {
      map_intrinsics[indexCam] = intrinsic_it.second->getParams();

      double * parameter_block = &map_intrinsics[indexCam][0];
      problem.AddParameterBlock(parameter_block, map_intrinsics[indexCam].size());
      if (options.intrinsics_opt == Intrinsic_Parameter_Type::NONE)
      {
        // set the whole parameter block as constant for best performance
        problem.SetParameterBlockConstant(parameter_block);
      }
      else
      {
        const std::vector<int> vec_constant_intrinsic =
          intrinsic_it.second->subsetParameterization(options.intrinsics_opt);
        if (!vec_constant_intrinsic.empty())
        {
          ceres::SubsetParameterization *subset_parameterization =
            new ceres::SubsetParameterization(
              map_intrinsics[indexCam].size(), vec_constant_intrinsic);
          problem.SetParameterization(parameter_block, subset_parameterization);
        }
      }
    }
    else
    {
      std::cerr << "Unsupported camera type." << std::endl;
    }
  }

  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction =
    ceres_options.bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (auto & structure_landmark_it : sfm_data.structure)
  {
    const Observations & obs = structure_landmark_it.second.obs;

    for (const auto & obs_it : obs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(obs_it.first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), obs_it.second.x);

      if (cost_function)
        problem.AddResidualBlock(cost_function,
          p_LossFunction,
          &map_intrinsics[view->id_intrinsic][0],
          &map_poses[view->id_pose][0],
          structure_landmark_it.second.X.data());
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(structure_landmark_it.second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (auto & gcp_landmark_it : sfm_data.control_points)
    {
      const Observations & obs = gcp_landmark_it.second.obs;

      for (const auto & obs_it : obs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(obs_it.first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            obs_it.second.x,
            options.control_point_opt.weight);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            nullptr,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            gcp_landmark_it.second.X.data());
      }
      if (obs.empty())
      {
        std::cerr
          << "Cannot use this GCP id: " << gcp_landmark_it.first
          << ". There is not linked image observation." << std::endl;
      }
      else
      {
        // Set the 3D point as FIXED (it's a valid GCP)
        problem.SetParameterBlockConstant(gcp_landmark_it.second.X.data());
      }
    }
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
  ceres::Solver::Options ceres_config_options;
  ceres_config_options.max_num_iterations = 500;
  ceres_config_options.preconditioner_type = ceres_options.preconditioner_type_;
  ceres_config_options.linear_solver_type = ceres_options.linear_solver_type_;
  ceres_config_options.sparse_linear_algebra_library_type = ceres_options.sparse_linear_algebra_library_type_;
  ceres_config_options.minimizer_progress_to_stdout = ceres_options.bVerbose_;
  ceres_config_options.logging_type = ceres::SILENT;
  ceres_config_options.num_threads = ceres_options.nb_threads_;
  ceres_config_options.num_linear_solver_threads = ceres_options.nb_threads_;
  ceres_config_options.parameter_tolerance = ceres_options.parameter_tolerance_;

  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(ceres_config_options, &problem, &summary);
  if (ceres_options.bCeres_summary_)
    std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (!summary.IsSolutionUsable())
  {
    if (ceres_options.bVerbose_)
      std::cout << "Bundle Adjustment failed." << std::endl;
    return false;
  }
  else // Solution is usable
  {
    if (ceres_options.bVerbose_)
    {
      // Display statistics about the minimization
      std::cout << std::endl
        << "Bundle Adjustment statistics (approximated RMSE):\n"
        << " #views: " << sfm_data.views.size() << "\n"
        << " #poses: " << sfm_data.poses.size() << "\n"
        << " #intrinsics: " << sfm_data.intrinsics.size() << "\n"
        << " #tracks: " << sfm_data.structure.size() << "\n"
        << " #residuals: " << summary.num_residuals << "\n"
        << " Initial RMSE: " << std::sqrt( summary.initial_cost / summary.num_residuals) << "\n"
        << " Final RMSE: " << std::sqrt( summary.final_cost / summary.num_residuals) << "\n"
        << " Time (s): " << summary.total_time_in_seconds << "\n"
        << std::endl;
    }

    // Update camera poses with refined data
    if (options.extrinsics_opt != Extrinsic_Parameter_Type::NONE)
    {
      for (auto & pose_it : sfm_data.poses)
      {
        const IndexT indexPose = pose_it.first;

        Mat3 R_refined;
        ceres::AngleAxisToRotationMatrix(&map_poses[indexPose][0], R_refined.data());
        Vec3 t_refined(map_poses[indexPose][3], map_poses[indexPose][4], map_poses[indexPose][5]);
        // Update the pose
        Pose3 & pose = pose_it.second;
        pose = Pose3(R_refined, -R_refined.transpose() * t_refined);
      }
    }

    // Update camera intrinsics with refined data
    if (options.intrinsics_opt != Intrinsic_Parameter_Type::NONE)
    {
      for (auto & intrinsic_it : sfm_data.intrinsics)
      {
        const IndexT indexCam = intrinsic_it.first;

        const std::vector<double> & vec_params = map_intrinsics[indexCam];
        intrinsic_it.second.get()->updateFromParams(vec_params);
      }
    }
    
    size_t nbOutliers_residualErr;
    size_t nbOutliers_angleErr;
    do
    {
      nbOutliers_residualErr = RemoveOutliers_PixelResidualError(sfm_data, 4.0, 2);
      nbOutliers_angleErr = RemoveOutliers_AngleError(sfm_data, 2.0);
    }while((nbOutliers_residualErr + nbOutliers_angleErr) > 0);
    
    return true;
  }
}




// ----------------------------------------------------
// Localization of the SfM reconstruction in another SfM reconstruction (and alignment)
// ----------------------------------------------------
int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "  Localization of the SfM reconstruction in another SfM reconstruction (and alignment):\n"
    << "-----------------------------------------------------------\n"
    << std::endl;

  CmdLine cmd;

  std::string sSfM_Data_Filename_A;
  std::string sSfM_Data_Filename_B;
  std::string sMatchesDir_A;
  std::string sMatchesDir_B;
  std::string sOutDir = "";
  int iIntrinsicsTransferMode = 0;
  double dMaxResidualError = std::numeric_limits<double>::infinity();
  double dPoseResidualMedianMultiply = 10;
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";

#ifdef OPENMVG_USE_OPENMP
  int iNumThreads = 0;
#endif

  cmd.add( make_option('i', sSfM_Data_Filename_A, "input_file_A") );
  cmd.add( make_option('j', sSfM_Data_Filename_B, "input_file_B") );
  cmd.add( make_option('a', sMatchesDir_A, "match_dir_A") );
  cmd.add( make_option('b', sMatchesDir_B, "match_dir_B") );
  cmd.add( make_option('o', sOutDir, "out_dir") );
  cmd.add( make_option('t', iIntrinsicsTransferMode, "intrinsics_mode") );
  
  cmd.add( make_option('r', dMaxResidualError, "residual_error"));
  cmd.add( make_option('d', dPoseResidualMedianMultiply, "pose_residual_error"));
  cmd.add( make_option('f', sIntrinsic_refinement_options, "refineIntrinsics") );
#ifdef OPENMVG_USE_OPENMP
  cmd.add( make_option('n', iNumThreads, "numThreads") );
#endif

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--input_file_A] path to base SfM_Data scene (origin)\n"
    << "[-j|--input_file_B] path to modifying SfM_Data scene\n"
    << "[-a|--match_dir_A] path to the directory containing the matches\n"
    << "  corresponding to the base SfM_Data scene (specified with -i)\n"
    << "[-b|--match_dir_B] path to the directory containing the matches\n"
    << "  corresponding to the second SfM_Data scene (specified with -j)\n"
    << "[-o|--out_dir] path where the output data will be stored\n"
    << "[-t|--intrinsics_mode] Which intrinsic the redifined reconstruction uses:\n"
      << "\t 0: As defined in its own reconstruction \n"
      << "\t 1: As defined in base reconstruction (default)\n"
      << "\t\t (only possible if a single intrinsic is provided in base reconstruction)\n"
    << "\n"
    << "(optional)\n"
    << "[-r|--residual_error] upper bound of the residual error tolerance\n"
    << "[-d|--pose_residual_error] multiplication factor of median error for outlier threshold in pose localization\n"
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
#ifdef OPENMVG_USE_OPENMP
    << "[-n|--numThreads] number of thread(s)\n"
#endif
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }



  // Load input SfM_Data scene
  // Set A
  SfM_Data sfm_data_A;
  if (!Load(sfm_data_A, sSfM_Data_Filename_A, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data_A file \""<< sSfM_Data_Filename_A << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if (sfm_data_A.GetPoses().empty() || sfm_data_A.GetLandmarks().empty())
  {
    std::cerr << std::endl
      << "The input SfM_Data_A file have not 3D content to match with." << std::endl;
    return EXIT_FAILURE;
  }
  
  // Set B
  SfM_Data sfm_data_B;
  if (!Load(sfm_data_B, sSfM_Data_Filename_B, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data_B file \""<< sSfM_Data_Filename_B << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }

  if (sfm_data_B.GetPoses().empty() || sfm_data_B.GetLandmarks().empty())
  {
    std::cerr << std::endl
      << "The input SfM_Data_B file have not 3D content to match with." << std::endl;
    return EXIT_FAILURE;
  }

  // ---------------
  // Initialization
  // ---------------

  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  // Set A
  const std::string sImage_describer_A = stlplus::create_filespec(sMatchesDir_A, "image_describer", "json");
  std::unique_ptr<Regions> regions_type_A = Init_region_type_from_file(sImage_describer_A);
  if (!regions_type_A)
  {
    std::cerr << "Invalid: "
      << sImage_describer_A << " regions type file." << std::endl;
    return EXIT_FAILURE;
  }
  // Set B
  const std::string sImage_describer_B = stlplus::create_filespec(sMatchesDir_B, "image_describer", "json");
  std::unique_ptr<Regions> regions_type_B = Init_region_type_from_file(sImage_describer_B);
  if (!regions_type_B)
  {
    std::cerr << "Invalid: "
      << sImage_describer_B << " regions type file." << std::endl;
    return EXIT_FAILURE;
  } 
  
  
  // Init the feature extractor that have been used for the reconstruction
  // Set A
  std::unique_ptr<Image_describer> image_describer_A;
  if (stlplus::is_file(sImage_describer_A))
  {
    // Dynamically load the image_describer from the file (will restore old used settings)
    std::ifstream stream(sImage_describer_A.c_str());
    if (!stream.is_open())
      return false;

    try
    {
      cereal::JSONInputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", image_describer_A));
    }
    catch (const cereal::Exception & e)
    {
      std::cerr << e.what() << std::endl
        << "Cannot dynamically allocate the Image_describer interface." << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    std::cerr << "Expected file image_describer.json cannot be opened." << std::endl;
    return EXIT_FAILURE;
  }
  // Set B
  std::unique_ptr<Image_describer> image_describer_B;
  if (stlplus::is_file(sImage_describer_B))
  {
    // Dynamically load the image_describer from the file (will restore old used settings)
    std::ifstream stream(sImage_describer_B.c_str());
    if (!stream.is_open())
      return false;

    try
    {
      cereal::JSONInputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", image_describer_B));
    }
    catch (const cereal::Exception & e)
    {
      std::cerr << e.what() << std::endl
        << "Cannot dynamically allocate the Image_describer interface." << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    std::cerr << "Expected file image_describer.json cannot be opened." << std::endl;
    return EXIT_FAILURE;
  }  
  
  // Load the SfM_Data region's views
  std::shared_ptr<Regions_Provider> regions_provider_A = std::make_shared<Regions_Provider>();
  
  if (!regions_provider_A->load(sfm_data_A, sMatchesDir_A, regions_type_A)) {
    std::cerr << std::endl << "Invalid regions." << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);


  if (iIntrinsicsTransferMode < 0 ||
      iIntrinsicsTransferMode > 1 )  {
    std::cerr << "\n Invalid intrinsic transfer mode" << std::endl;
    return EXIT_FAILURE;
  }  

    // Check if there is only one intrinsics in base reconstruction
  if(iIntrinsicsTransferMode == 1 && (sfm_data_A.intrinsics.size() != 1 || sfm_data_B.intrinsics.size() != 1) )
  {
    // If not set transfer mode to estimate while localizing
    std::cout<<"Both reconstructions do not have a single intrinsics!"<<std::endl
    << "Using original intrinsics of second reconstruction."<<std::endl;
    
    iIntrinsicsTransferMode = 0;
  }

  
  // Initialization of retrieval database
  sfm::SfM_Localization_Single_3DTrackObservation_Database localizer;
  if (!localizer.Init(sfm_data_A, *regions_provider_A.get()))
  {
    std::cerr << "Cannot initialize the SfM localizer" << std::endl;
  }
  // Since we have copied interesting data, release some memory
  regions_provider_A.reset();
  
  

  // Estimated poses of views in the reconstruction B
  Hash_Map<IndexT, Pose3> poses_estimated;

#ifdef OPENMVG_USE_OPENMP
  const unsigned int nb_max_thread = (iNumThreads == 0) ? 0 : omp_get_max_threads();
  omp_set_num_threads(nb_max_thread);
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < static_cast<int>(sfm_data_B.poses.size()); ++i)
  {
    // For each pose in reconstruction B try to localize
    Poses::const_iterator pose_it = sfm_data_B.poses.begin();
    std::advance(pose_it, i);
    
    const IndexT indexPose_B = pose_it->first; // Also view id
    const View * view_B = sfm_data_B.views.at(indexPose_B).get();
    
    std::unique_ptr<Regions> query_regions(regions_type_B->EmptyClone());
    {
      const std::string &sView_filename = view_B->s_Img_path;

      // Read features 
      const std::string
        sFeat = stlplus::create_filespec(sMatchesDir_B, stlplus::basename_part(sView_filename.c_str()), "feat"),
        sDesc = stlplus::create_filespec(sMatchesDir_B, stlplus::basename_part(sView_filename.c_str()), "desc");

      // Skip views whos features we cant read
      if (!stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc))
      {
        continue;
      }
      else // load already existing regions
      {
        query_regions->Load(sFeat,sDesc);
      }
    }
    
    std::shared_ptr<cameras::IntrinsicBase> optional_intrinsic(nullptr);
    
    switch(iIntrinsicsTransferMode)
    {
      case 0:
        optional_intrinsic = sfm_data_B.intrinsics[view_B->id_intrinsic];
        break;
      case 1:
        optional_intrinsic = sfm_data_A.intrinsics[0];
        break;
    }
    
    Pose3 pose_estimated;
    sfm::Image_Localizer_Match_Data matching_data;
    matching_data.error_max = dMaxResidualError;
    
    const IndexT &view_B_width = view_B->ui_width;
    const IndexT &view_B_height = view_B->ui_height;
    // Try to localize the view
    if (!localizer.Localize(
      Pair(view_B_width, view_B_height),
      optional_intrinsic.get(),
      *(query_regions.get()),
      pose_estimated,
      &matching_data))
    {
      std::cerr << "Cannot locate the image " << view_B->s_Img_path << std::endl;
    }
    else
    {
      // A valid pose has been found (try to refine it):
      if (!sfm::SfM_Localizer::RefinePose
      (
        optional_intrinsic.get(),
        pose_estimated, matching_data,
        true, false
      ))
      {
        std::cerr << "Refining pose for image " << view_B->s_Img_path << " failed." << std::endl;
      }
      
      // Save the new pose
#ifdef OPENMVG_USE_OPENMP
      #pragma omp critical
#endif
      poses_estimated[view_B->id_view] = pose_estimated;
    }
  }
  
  std::cout << " Total poses found : " << poses_estimated.size() << "/" << sfm_data_B.poses.size() << endl;
  
  // Check for outliers
  if (sfm_data_B.GetPoses().size() > 3)
  {
    double dEstimatedPose_Outlier_Threshold;
  
    // Compute the registration:
    Similarity3 sim;
    std::vector<Vec3> X_B, X_B_A;
    // Add original and new estimated poses
    for (const auto & view_it : poses_estimated)
    {
      const IndexT view_id = view_it.first;    
      X_B.push_back( sfm_data_B.GetPoses().at(view_id).center() );
      X_B_A.push_back(view_it.second.center() );
    }
    
    Mat X_B_Mat = Eigen::Map<Mat>(X_B[0].data(),3, X_B.size());
    Mat X_B_A_Mat = Eigen::Map<Mat>(X_B_A[0].data(),3, X_B_A.size());
    kernel::Similarity3_Kernel kernel_noisy(X_B_Mat, X_B_A_Mat);
    
    // Estimate similarty transformation between centers of cameras
    double lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel_noisy, &sim);
    
    std::cout<<"Noisy similarity transform median: "<<lmeds_median<<"\n";
    if (lmeds_median == std::numeric_limits<double>::max())
    {
      std::cerr << "Second reconstruction cannot be accurately localized!"<< std::endl
      <<"\t - Inaccurate localization of views"<< std::endl;
      return EXIT_FAILURE;
    }
    
    // Compute the median residual error once the registration is applied
    for (Vec3 & pos : X_B) // Transform B poses for residual computation
    {
      pos = sim(pos);
    }
    Vec residual = (Eigen::Map<Mat3X>(X_B[0].data(), 3, X_B.size()) - Eigen::Map<Mat3X>(X_B_A[0].data(), 3, X_B_A.size())).colwise().norm();
    std::sort(residual.data(), residual.data() + residual.size());
    dEstimatedPose_Outlier_Threshold = residual(residual.size()/2) * dPoseResidualMedianMultiply;
          
    // Remove poses with large residual error (probable outliers)  
    Hash_Map<IndexT, Pose3>::iterator pose_it = poses_estimated.begin();
    while(pose_it != poses_estimated.end())
    {
      IndexT view_id = pose_it->first; 
      double residual = (sim(sfm_data_B.GetPoses().at(view_id).center()) - pose_it->second.center()).norm();
      if(residual > dEstimatedPose_Outlier_Threshold)
      {
        pose_it = poses_estimated.erase(pose_it);
      }
      else
      {
        ++pose_it;
      }
    }

    std::cout << " Total poses after outlier removal : " << poses_estimated.size() << std::endl;

    // Recompute similarity transformation
    X_B.clear();
    X_B_A.clear();
    for (const auto & view_it : poses_estimated)
    {
      const IndexT view_id = view_it.first;    
      X_B.push_back( sfm_data_B.GetPoses().at(view_id).center() );
      X_B_A.push_back(view_it.second.center() );
    }
    
    X_B_Mat = Eigen::Map<Mat>(X_B[0].data(),3, X_B.size());
    X_B_A_Mat = Eigen::Map<Mat>(X_B_A[0].data(),3, X_B_A.size());
    kernel::Similarity3_Kernel kernel(X_B_Mat, X_B_A_Mat);
    lmeds_median = openMVG::robust::LeastMedianOfSquares(kernel, &sim);
    
    if (lmeds_median == std::numeric_limits<double>::max())
    {
      std::cerr << "Second reconstruction cannot be accurately localized!"<< std::endl
      <<"\t - Inaccurate localization of views"<< std::endl;
      return EXIT_FAILURE;
    }
    
    std::cout<<"Similarity transform median: "<<lmeds_median<<"\n";
    
    // Go through poses
    // the ones with good localization use it otherwise approximate with similarity
    
    for (auto & pose_it : sfm_data_B.poses)
    {
      IndexT indexPose = pose_it.first;
      if(poses_estimated.count(indexPose) != 0)
      {
        sfm_data_B.poses[indexPose] = poses_estimated[indexPose];
      }
      else
      {
        sfm_data_B.poses[indexPose] = sim(sfm_data_B.poses[indexPose]);
      }
    }
      
    for (auto & landmark_it : sfm_data_B.structure)
    {
      landmark_it.second.X = sim(landmark_it.second.X);
    }
    for (auto & landmark_it : sfm_data_B.control_points)
    {
      landmark_it.second.X = sim(landmark_it.second.X);
    }
  }
  else
  {
      std::cerr << "Second reconstruction cannot be accurately localized!"<< std::endl
      <<"\t - Insufficient number of localized views"<< std::endl;
      
      return EXIT_FAILURE;
  }
  
  // If we use one intrinsics for all (from base reconstruction) we save it
  if(iIntrinsicsTransferMode == 1)
  {
    sfm_data_B.intrinsics[0] = sfm_data_A.intrinsics[0];
  }
  
  GroupSharedIntrinsics(sfm_data_B);
  
  std::cout << " Total poses after outlier removal : " << poses_estimated.size() << "/" << sfm_data_B.poses.size() << endl;
  
  
  // Perform BA to compensate for new changes in poses and intrinsics
  
  // If we use camera params of first reconstruction the intrinsics are fixed
  if(iIntrinsicsTransferMode == 1)
  {
    sIntrinsic_refinement_options = "NONE";
  }
  const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options = cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);
  
  // Adjust the cameras that were not sucessfully matched and the structure
  AdjustReconstructionWithSomeFixedCameras(sfm_data_B, intrinsic_refinement_options, poses_estimated);
  
  
  // Export the found camera position in a ply.
  //const std::string out_file_name = stlplus::create_filespec(sOutDir, "found_pose_centers", "ply");
  //plyHelper::exportToPly(vec_found_poses, out_file_name);

  // Export found camera poses along with original reconstruction in a new sfm_data file
  if (!Save(
    sfm_data_B,
    stlplus::create_filespec( sOutDir, "sfm_data.bin" ).c_str(),
    ESfM_Data(ALL)))
  {
    return EXIT_FAILURE;
  }
  // export also as ply
  if (!Save(
    sfm_data_B,
    stlplus::create_filespec( sOutDir, "cloud_and_points_aligned.ply" ).c_str(),
    ESfM_Data(ALL)))
  {
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
