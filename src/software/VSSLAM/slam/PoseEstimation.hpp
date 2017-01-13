
#pragma once

#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/solver_homography_kernel.hpp"

#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"


namespace openMVG  {
namespace VSSLAM  {

using namespace openMVG::robust;

static const size_t ACRANSAC_ITER = 64;

static double computeHomographyScore
(
    const Mat & model,
    const Mat & x1,
    const Mat & x2,
    const double gamma
)
{
  Mat model_inv = model.inverse();

  double score = 0.0;
  /*
  Mat d2_1in2 = (x2 - ((model * x1.colwise().homogeneous()).colwise().hnormalized())).colwise().squaredNorm();
  Mat d2_2in1 = (x1 - ((model_inv * x2.colwise().homogeneous()).colwise().hnormalized())).colwise().squaredNorm();

  for (size_t pt_i = 0; pt_i < x1.cols(); ++pt_i)
  {
    const double d2_1_2 = d2_1in2(pt_i);
    const double d2_2_1 = d2_2in1(pt_i);
    if (d2_1_2 < gamma)
      score += d2_1_2;
    if (d2_2_1 < gamma)
      score += d2_2_1;
  }*/

  for (size_t pt_i = 0; pt_i < x1.cols(); ++pt_i)
  {
    const Vec2 & pt_1 = x1.col(pt_i);
    const Vec2 & pt_2 = x2.col(pt_i);
    // x2 = H * x1
    Vec2 proj_x = Vec3(model * pt_1.homogeneous()).hnormalized();
    double d2 = (pt_2 - proj_x).squaredNorm();

    if (d2 < gamma)
      score += (gamma - d2);

    // x1 = H^-1 * x2 hnormalized()
    proj_x = Vec3(model_inv * pt_2.homogeneous()).hnormalized();
    d2 = (pt_1 - proj_x).squaredNorm();

    if (d2 < gamma)
      score += (gamma - d2);

  }

  return score;
}

static double computeEpipolarScore
(
    const Mat & model,
    const Mat & x1,
    const Mat & x2,
    const double gamma
)
{
  Mat model_inv = model.inverse();

  double score = 0.0;

  /*
  Mat FFX_1 = model * x1.colwise().homogeneous();
  Mat d2_1in2 = ((FFX_1.cwiseProduct(x2.colwise().homogeneous())).colwise().sum().array().pow(2.0f)) / (FFX_1.topRows(2).colwise().squaredNorm()).array();
  Mat FFX_2 = model_inv * x2.colwise().homogeneous();
  Mat d2_2in1 = ((FFX_2.cwiseProduct(x1.colwise().homogeneous())).colwise().sum().array().pow(2.0f)) / (FFX_2.topRows(2).colwise().squaredNorm()).array();

  for (size_t pt_i = 0; pt_i < x1.cols(); ++pt_i)
  {
    const double d2_1_2 = d2_1in2(pt_i);
    const double d2_2_1 = d2_2in1(pt_i);
    if (d2_1_2 < gamma)
      score += d2_1_2;
    if (d2_2_1 < gamma)
      score += d2_2_1;
  }*/

  for (size_t pt_i = 0; pt_i < x1.cols(); ++pt_i)
  {
    const Vec3 pt_2 = x2.col(pt_i).homogeneous();
    const Vec3 pt_1 = x1.col(pt_i).homogeneous();

    // x2T * F * x1 = 0
    Vec3 F_x = model * pt_1;
    double dF_x = F_x.dot(pt_2);
    double d2 = (dF_x * dF_x) /  F_x.head<2>().squaredNorm();

    if (d2 < gamma)
      score += (gamma - d2);

    // x1T * F^-1 * x2 = 0
    F_x = model_inv * pt_2;
    dF_x = F_x.dot(pt_1);
    d2 = (dF_x * dF_x) /  F_x.head<2>().squaredNorm();

    if (d2 < gamma)
      score += (gamma - d2);

  }

  return score;
}


  static bool computeH
  (
    const Mat & x1,
    const Mat & x2,
    const unsigned int w1,
    const unsigned int h1,
    const unsigned int w2,
    const unsigned int h2,
    Mat3 & H_21,  //x2 = H21 * x1
    double & threshold
  )
  {
    //-- Homography robust estimation
    using KernelType =
      ACKernelAdaptor<
        openMVG::homography::kernel::FourPointSolver,
        openMVG::homography::kernel::AsymmetricError,
        UnnormalizerI,
        Mat3>;

    KernelType kernel(
      x1, w1, h1,
      x2, w2, h2,
      false); // configure as point to point error model.

    std::vector<size_t> vec_inliers;

    // ACRansacOut.first (threshold) ACRansacOut.second(NFA score)
    const std::pair<double,double> ACRansacOut = ACRANSAC(kernel, vec_inliers, ACRANSAC_ITER, &H_21, 4, false);

    std::cout<<"Computed H: "<<vec_inliers.size()<<" :: thresh: "<<ACRansacOut.first<<" :: NFA: "<<ACRansacOut.second<<"\n";

    if (vec_inliers.size() > KernelType::MINIMUM_SAMPLES *2.5)
    {
      threshold = ACRansacOut.first;
      return true;
    }
    return false;

    // Check the homography support some point to be considered as valid
    //return vec_inliers.size() > KernelType::MINIMUM_SAMPLES *2.5;
  }

  static bool computeH
  (
    const Mat & x1,
    const Mat & x2,
    const unsigned int w,
    const unsigned int h,
    Mat3 & H,  //x2 = H21 * x1
    double & threshold
  )
  {
    return computeH(x1, x2, w, h, w, h, H, threshold);
  }


  static bool computeE
  (
    const Mat3 & K1,
    const Mat3 & K2,
    const Mat & x1,
    const Mat & x2,
    const unsigned int w1,
    const unsigned int h1,
    const unsigned int w2,
    const unsigned int h2,
    Mat3 & pE,
    double & threshold
  )
  {
    // Use the 5 point solver to estimate E
    using SolverType = openMVG::essential::kernel::FivePointKernel;
    // Define the AContrario adaptor
    using KernelType =
      ACKernelAdaptorEssential<
        SolverType,
        openMVG::fundamental::kernel::EpipolarDistanceError,
        Mat3>;

    KernelType kernel(
      x1, w1, h1,
      x2, w2, w2,
      K1, K2);

    // Robustly estimation of the Essential matrix and it's precision
    std::vector<size_t> vec_inliers;
    // ACRansacOut.first (threshold) ACRansacOut.second(NFA score)
    const std::pair<double,double> ACRansacOut = ACRANSAC(kernel, vec_inliers, ACRANSAC_ITER, &pE, 4, false);

    std::cout<<"Computed E: "<<vec_inliers.size()<<" :: thresh: "<<ACRansacOut.first<<" :: NFA: "<<ACRansacOut.second<<"\n";

    if (vec_inliers.size() > SolverType::MINIMUM_SAMPLES *2.5)
    {
      threshold = ACRansacOut.first;
      return true;
    }
    return false;

    //return vec_inliers.size() > 2.5 * SolverType::MINIMUM_SAMPLES;
  }

  static bool computeE
  (
    const Mat3 & K,
    const Mat & x1,
    const Mat & x2,
    const unsigned int w,
    const unsigned int h,
    Mat3 & pE,
    double & threshold
  )
  {
    return computeE(K, K, x1, x2, w, h, w, h, pE, threshold);
  }


}
}
