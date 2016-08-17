// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_uncertainty.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres_camera_functor.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/types.hpp"

#include <cereal/cereal.hpp> // Serialization

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SVD>

namespace openMVG {
namespace sfm {

using namespace openMVG::cameras;
using namespace openMVG::geometry;



Eigen::MatrixXd pseudoInverse(const Eigen::MatrixXd &a, double epsilon = std::numeric_limits<double>::epsilon())
{
	Eigen::JacobiSVD< Eigen::MatrixXd > svd(a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().transpose();
}
Eigen::SparseMatrix<double, Eigen::RowMajor> pseudoInverse(const Eigen::SparseMatrix<double, Eigen::RowMajor> &a, double epsilon = std::numeric_limits<double>::epsilon()){
  Eigen::MatrixXd dense_a = Eigen::MatrixXd(a);

	Eigen::JacobiSVD< Eigen::MatrixXd > svd(dense_a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
  double tolerance = epsilon * std::max(dense_a.cols(), dense_a.rows()) *svd.singularValues().array().abs()(0);
  Eigen::MatrixXd dense_a_inv = svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().transpose();

  Eigen::SparseMatrix<double, Eigen::RowMajor> sparse_a_inv =  dense_a_inv.sparseView();
  sparse_a_inv.makeCompressed();
  return sparse_a_inv;
}


bool EstimateUncertainty
(
  SfM_Data & sfm_data,     // the SfM scene to refine
  Bundle_Adjustment_Ceres &bundle_adjustment_obj,
  const Optimize_Options &options,
  const bool evaluateLandmarks
)
{
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
 for (Poses::const_iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose)
  {
    const IndexT indexPose = itPose->first;

    const Pose3 & pose = itPose->second;
    const Mat3 R = pose.rotation();
    const Vec3 t = pose.translation();

    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // angleAxis + translation
    map_poses[indexPose] = {angleAxis[0], angleAxis[1], angleAxis[2], t(0), t(1), t(2)};

    double * parameter_block = &map_poses[indexPose][0];
    problem.AddParameterBlock(parameter_block, 6);
    if (options.extrinsics_opt == Extrinsic_Parameter_Type::NONE)
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
  for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
    itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
  {
    const IndexT indexCam = itIntrinsic->first;

    if (isValid(itIntrinsic->second->getType()))
    {
      map_intrinsics[indexCam] = itIntrinsic->second->getParams();

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
          itIntrinsic->second->subsetParameterization(options.intrinsics_opt);
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
    bundle_adjustment_obj.ceres_options().bUse_loss_function_ ?
      new ceres::HuberLoss(Square(4.0))
      : nullptr;

  // For all visibility add reprojections errors:
  for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks!= sfm_data.structure.end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;

    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      // Build the residual block corresponding to the track observation:
      const View * view = sfm_data.views.at(itObs->first).get();

      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.
      ceres::CostFunction* cost_function =
        IntrinsicsToCostFunction(sfm_data.intrinsics[view->id_intrinsic].get(), itObs->second.x);

      if (cost_function)
        problem.AddResidualBlock(cost_function,
          p_LossFunction,
          &map_intrinsics[view->id_intrinsic][0],
          &map_poses[view->id_pose][0],
          iterTracks->second.X.data());
    }
    if (options.structure_opt == Structure_Parameter_Type::NONE)
      problem.SetParameterBlockConstant(iterTracks->second.X.data());
  }

  if (options.control_point_opt.bUse_control_points)
  {
    // Use Ground Control Point:
    // - fixed 3D points with weighted observations
    for (Landmarks::iterator iterGCPTracks = sfm_data.control_points.begin();
      iterGCPTracks!= sfm_data.control_points.end(); ++iterGCPTracks)
    {
      const Observations & obs = iterGCPTracks->second.obs;

      for (Observations::const_iterator itObs = obs.begin();
        itObs != obs.end(); ++itObs)
      {
        // Build the residual block corresponding to the track observation:
        const View * view = sfm_data.views.at(itObs->first).get();

        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.
        ceres::CostFunction* cost_function =
          IntrinsicsToCostFunction(
            sfm_data.intrinsics[view->id_intrinsic].get(),
            itObs->second.x,
            options.control_point_opt.weight);

        if (cost_function)
          problem.AddResidualBlock(cost_function,
            nullptr,
            &map_intrinsics[view->id_intrinsic][0],
            &map_poses[view->id_pose][0],
            iterGCPTracks->second.X.data());
      }
      // Set the 3D point as FIXED (it's a GCP)
      problem.SetParameterBlockConstant(iterGCPTracks->second.X.data());
    }
  }


  // Confgure problem evaluator
  ceres::Problem::EvaluateOptions evaluate_options_;
  evaluate_options_.num_threads = bundle_adjustment_obj.ceres_options().nb_threads_;
  ceres::CRSMatrix jacobian;
  double cost;
  // Evaluate problem
  problem.Evaluate(evaluate_options_, &cost, NULL, NULL, &jacobian);

  // -----------------------------------------------
  // Convert Ceres Jacobian to Sparse Eigen Matrix
  // -----------------------------------------------
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor> EigenSparseMatrix;
  // Convert the matrix to column major order
  EigenSparseMatrix sparse_jacobian =
      Eigen::MappedSparseMatrix<double, Eigen::RowMajor>(
          jacobian.num_rows, jacobian.num_cols,
          static_cast<int>(jacobian.values.size()),
          jacobian.rows.data(), jacobian.cols.data(), jacobian.values.data());
          
  // -----------------------------------------------
  // Compute size of parameters
  // -----------------------------------------------
  // Compute the size of different parameter blocks
  const int total_observations = sparse_jacobian.rows();
  const int total_cam_ext_param = 6*sfm_data.poses.size();
  const int total_landmark_param = 3*sfm_data.structure.size();
  const int total_control_param = 3*sfm_data.control_points.size();
  int total_intrinsic_param = 0;
  // Save the index of the start of next intrinsic block of params
  Hash_Map<IndexT, int> start_intrinsic_param_per_row;
  for (Intrinsics::const_iterator itIntrinsic = sfm_data.intrinsics.begin();
    itIntrinsic != sfm_data.intrinsics.end(); ++itIntrinsic)
  {
    const int n_intrinsic_param = map_intrinsics[itIntrinsic->first].size();
    // Save the start index of the intrinsics
    start_intrinsic_param_per_row[itIntrinsic->first] = total_intrinsic_param;
    // Compute total number of intrinsics
    total_intrinsic_param += n_intrinsic_param;
  }
    
  // -----------------------------------------------
  // Extract Jacobian matrix relating cameras and landmarks
  // -----------------------------------------------
  EigenSparseMatrix J_A = sparse_jacobian.block(0,0,total_observations,total_cam_ext_param + total_intrinsic_param);
  EigenSparseMatrix J_B = sparse_jacobian.block(0,total_cam_ext_param + total_intrinsic_param,total_observations,total_landmark_param + total_control_param);
  
  
  // -----------------------------------------------
  // Compute E_x_inv - Covariance matrix of detected features
  // -----------------------------------------------
  if (bundle_adjustment_obj.ceres_options().bVerbose_)
  {
    std::cout<<"Computing E_x_inv\n";
  }

  EigenSparseMatrix E_x_inv(total_observations,total_observations);
  
  {
  // Assume all features with same covariance -> 2px
  Eigen::Matrix2d Ex;
  Ex << 0.5,0,0,0.5;
  IndexT obsID=0;
  for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
    iterTracks!= sfm_data.structure.end(); ++iterTracks)
  {
    const Observations & obs = iterTracks->second.obs;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      
      // Check if covariance of detected features is present
      if(sfm_data.uncertainty_observations.find(itObs->second.id_feat)== sfm_data.uncertainty_observations.end()){ 
        for(int r=0;r<2;r++){
          for(int c=0;c<2;c++){
            E_x_inv.insert(obsID*2+r,obsID*2+c) = Ex(r,c);
          }
        }
      }
      else{
        for(int r=0;r<2;r++){
          for(int c=0;c<2;c++){
            E_x_inv.insert(obsID*2+r,obsID*2+c) = sfm_data.uncertainty_observations[itObs->second.id_feat].covariance(r,c);            

          }
        }
      }
      obsID++;
    }
  }
  }
  
  // -----------------------------------------------
  // Compute U
  // ----------------------------------------------- 
  if (bundle_adjustment_obj.ceres_options().bVerbose_)
  { 
    std::cout<<"Computing U\n";  
  }
  EigenSparseMatrix U = (J_A.transpose() * E_x_inv * J_A).pruned(0.0);

  // -----------------------------------------------
  // Compute V_inverse
  // ----------------------------------------------- 
  if (bundle_adjustment_obj.ceres_options().bVerbose_)
  { 
    std::cout<<"Computing V_inverse\n";
  }
  EigenSparseMatrix V_inv = (J_B.transpose() * E_x_inv * J_B).pruned(0.0);
  
  // Invert diagonal blocks
  {
  Eigen::Matrix3d diag;
  for(int t_id=0;t_id<V_inv.rows();t_id+=3){
      diag = V_inv.block(t_id,t_id,3,3);
      diag = diag.inverse().eval();      
      for(int r=0;r<3;r++){
        for(int c=0;c<3;c++){
          V_inv.coeffRef(t_id+r,t_id+c) = diag(r,c);
        }
      }
  }
  }

  // -----------------------------------------------
  // Compute W
  // ----------------------------------------------- 
  if (bundle_adjustment_obj.ceres_options().bVerbose_)
  { 
    std::cout<<"Computing W\n";
  }  
  EigenSparseMatrix W = (J_A.transpose() * E_x_inv * J_B).pruned(0.0);

  // -----------------------------------------------
  // Compute Y and W * V_inverse * W'
  // ----------------------------------------------- 
  if (bundle_adjustment_obj.ceres_options().bVerbose_)
  { 
    std::cout<<"Computing Y\n";
  }
  EigenSparseMatrix Y = (W * V_inv);

  if (bundle_adjustment_obj.ceres_options().bVerbose_)
  {
    std::cout<<"Computing W * V_inv * W'\n";
  }
  EigenSparseMatrix WVW = (Y * W.transpose()).pruned(0.0);
  
  // -----------------------------------------------
  // Compute E_A
  // ----------------------------------------------- 
  if (bundle_adjustment_obj.ceres_options().bVerbose_)
  { 
    std::cout<<"Computing E_Poses\n";
  }  
  
  EigenSparseMatrix E_C = pseudoInverse(EigenSparseMatrix(U-WVW));
  int pose_seq=0;
  for (Poses::iterator itPose = sfm_data.poses.begin(); itPose != sfm_data.poses.end(); ++itPose,++pose_seq)
  {
    const IndexT indexPose = itPose->first;
    PoseUncertainty uncert_pose(indexPose);
    
    Eigen::Matrix<double, 6, 6> PP = E_C.block(6*pose_seq,6*pose_seq,6,6);
    uncert_pose.covariance = E_C.block(6*pose_seq,6*pose_seq,6,6);    
    sfm_data.uncertainty_poses.insert(std::make_pair(indexPose, uncert_pose));
  }


  // -----------------------------------------------
  // Compute E_P
  // -----------------------------------------------
  if(evaluateLandmarks){
    if (bundle_adjustment_obj.ceres_options().bVerbose_)
    { 
      std::cout<<"Computing E_Points\n";
    }

    EigenSparseMatrix W_column,WEW;
    Eigen::Matrix3d V_inv_i;
    // Loop through landmarks and compute its variance
    IndexT track_id=0;
    IndexT obs_id=0;
    for (Landmarks::iterator iterTracks = sfm_data.structure.begin();
      iterTracks!= sfm_data.structure.end(); ++iterTracks, ++track_id)
    {
      // Get i-th column of W (used for computing covariance
      W_column= W.middleCols(track_id*3,3);
      // i-th block of V
      V_inv_i = V_inv.block(track_id*3,track_id*3,3,3);
      // Compute E_P
      WEW = W_column.transpose() * E_C * W_column;

      LandmarkUncertainty uncert_point(iterTracks->first);
      uncert_point.covariance = (V_inv_i * WEW * V_inv_i) + V_inv_i;
      sfm_data.uncertainty_structure.insert(std::make_pair(iterTracks->first,uncert_point));
    }
  }
  return true;
}


/// Export uncertainty of points as double vector
void EstimateQualityOfStructure(const SfM_Data & sfm_data, std::vector<double> & vec_structureUncertainty)
{
  IndexT cpt = 0;
  for (Landmarks::const_iterator it = sfm_data.GetLandmarks().begin();
    it != sfm_data.GetLandmarks().end(); ++it, ++cpt)
  {
    double quality=0.0;
    
    // Compute the confidence (95%) ellipsoid and compute its volume
    if((sfm_data.uncertainty_structure.find(it->first)!=sfm_data.uncertainty_structure.end())){
        Eigen::Matrix3d mat = (sfm_data.GetLandmarksUncertainty().at(it->first).covariance);
Eigen::Vector3d eigen_values = mat.eigenvalues().real();
      // Parameters of ellipsoid
      double a_ellipsoid = sqrt(eigen_values(0)*7.815)*2;
      double b_ellipsoid = sqrt(eigen_values(1)*7.815)*2;
      double c_ellipsoid = sqrt(eigen_values(2)*7.815)*2;
      // Volume
      quality = (1.33333f * M_PI * a_ellipsoid * b_ellipsoid  * c_ellipsoid);
    }
    
    /*
    // Trace of covariance matrix
    double quality = (it->second.covariance).trace();
    //double qualityB = sfm_data.structure.at(cpt).covariance.trace();
    std::cout<<"TT: "<<quality<<" :: "<<"\n";
        //* (it->second.meanReprojError);
    */

    double meanReprojError = 0.0;
    const Observations & obs = it->second.obs;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      const Observation obs = itObs->second;
      const IndexT view_id = itObs->first;
      const View * view = sfm_data.GetViews().at(view_id).get();
      const IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      meanReprojError += cam->residual(pose,it->second.X,cam->get_ud_pixel(itObs->second.x)).squaredNorm();
    }
    meanReprojError /=obs.size();
  
    quality = quality * meanReprojError;

    vec_structureUncertainty.push_back(quality);
  }
}

} // namespace sfm
} // namespace openMVG

