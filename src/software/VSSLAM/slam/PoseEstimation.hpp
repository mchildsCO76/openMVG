
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

  static bool computeH
  (
    const Mat & x1,
    const Mat & x2,
    const std::pair<size_t, size_t> & size_ima,
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
      x1, size_ima.first, size_ima.second,
      x2, size_ima.first, size_ima.second,
      false); // configure as point to point error model.

    std::vector<size_t> vec_inliers;

    // ACRansacOut.first (threshold) ACRansacOut.second(NFA score)
    const std::pair<double,double> ACRansacOut = ACRANSAC(kernel, vec_inliers, 1024, &H_21, std::numeric_limits<double>::infinity(),
      true);

    std::cout<<"Computed H: "<<vec_inliers.size()<<" :: thresh: "<<ACRansacOut.first<<" :: NFA: "<<ACRansacOut.second<<"\n";

    /*
    if (vec_inliers.size() > KernelType::MINIMUM_SAMPLES *2.5)
    {
      std::vector<float> vec_errors (vec_inliers.size());
      Mat3 H_12 = H_21.inverse();

      for (size_t sample = 0; sample < vec_inliers.size(); ++sample)
      {
        float d2_1in2 = (x2.col(vec_inliers[sample]) - Vec3( H * x1.col(vec_inliers[sample]).homogeneous()).hnormalized() ).norm();
        float d2_2in1 = (x1.col(vec_inliers[sample]) - Vec3( H_inv * x2.col(vec_inliers[sample]).homogeneous()).hnormalized() ).norm();

        std::cout<<"ERROR "<<sample<<": "
            <<x1.col(vec_inliers[sample]).x()<<", "
            <<x1.col(vec_inliers[sample]).y()
            <<" : "
            <<x2.col(vec_inliers[sample]).x()<<", "
            <<x2.col(vec_inliers[sample]).y()
            <<" :: 1_2: "
            <<d2_1in2<<" :: 2_1: "
            <<d2_2in1<<"\n";

        //vec_errors[sample] = openMVG::homography::kernel::AsymmetricError::Error(H, x1.col(vec_inliers[sample]), x2.col(vec_inliers[sample]));
        //std::cout<<vec_errors[sample]<<" :: "<<err<<"\n";
      }
      score = 0.0;
      return true;
    }
    else
    {
      score = -1;
      return false;
    }*/

    if (vec_inliers.size() > KernelType::MINIMUM_SAMPLES *2.5)
    {
      threshold = ACRansacOut.first;
      return true;
    }
    return false;

    // Check the homography support some point to be considered as valid
    //return vec_inliers.size() > KernelType::MINIMUM_SAMPLES *2.5;
  }


  static bool computeE
  (
    const Mat3 & K,
    const Mat & x1,
    const Mat & x2,
    const std::pair<size_t, size_t> & size_ima,
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

    KernelType kernel(x1, size_ima.first, size_ima.second,
      x2, size_ima.first, size_ima.second, K, K);

    std::vector<size_t> vec_inliers;

    double precision = std::numeric_limits<double>::infinity();

    // Robustly estimation of the Essential matrix and it's precision
    // ACRansacOut.first (threshold) ACRansacOut.second(NFA score)
    const std::pair<double,double> ACRansacOut = ACRANSAC(kernel, vec_inliers, ACRANSAC_ITER, &pE, precision, false);

    if (vec_inliers.size() > SolverType::MINIMUM_SAMPLES *2.5)
    {
      threshold = ACRansacOut.first;
      return true;
    }
    return false;

    //return vec_inliers.size() > 2.5 * SolverType::MINIMUM_SAMPLES;
  }



}
}
