
// Copyright (c) 2016 Klemen ISTENIC.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once


#define _USE_MATH_DEFINES
#include <cmath>

#include <openMVG/types.hpp>

namespace openMVG {
namespace VSSLAM {

// Taken from Ceres AngleAxisToRotationMatrix (rotation.h line 288)
static void AngleAxisToRotationMatrix(const Vec3 & angle_axis, Mat3 & R) {
  const double theta2 = angle_axis.squaredNorm();
  if (theta2 > 0.0) {
    // We want to be careful to only evaluate the square root if the
    // norm of the angle_axis vector is greater than zero. Otherwise
    // we get a division by zero.
    double theta = sqrt(theta2);
    Vec3 w = angle_axis / theta;

    double costheta = cos(theta);
    double sintheta = sin(theta);

    R << costheta + w(0)*w(0) * (1.0 - costheta),
        w(2)*sintheta + w(0)*w(1) * (1.0 - costheta),
        -w(1)*sintheta + w(0)*w(2) * (1.0 - costheta),

        w(0)*w(1)*(1.0 - costheta) - w(2) * (sintheta),
        costheta + w(1) * w(1) * (1.0 - costheta),
        w(0)*sintheta + w(1) * w(2) * (1.0 - costheta),

        w(1)*sintheta + w(0) * w(2) * (1.0 - costheta),
        -w(0)*sintheta + w(1) * w(2) * (1.0 - costheta),
        costheta + w(2) * w(2) * (1.0 - costheta);


    /*R[0] =     costheta   + wx*wx*(1.0 -    costheta);
    R[1] =  wz*sintheta   + wx*wy*(1.0 -    costheta);
    R[2] = -wy*sintheta   + wx*wz*(1.0 -    costheta);

    R[3] =  wx*wy*(1.0 - costheta)     - wz*sintheta;
    R[4] =     costheta   + wy*wy*(1.0 -    costheta);
    R[5] =  wx*sintheta   + wy*wz*(1.0 -    costheta);

    R[6] =  wy*sintheta   + wx*wz*(1.0 -    costheta);
    R[7] = -wx*sintheta   + wy*wz*(1.0 -    costheta);
    R[8] =     costheta   + wz*wz*(1.0 -    costheta);*/
  } else {
    // At zero, we switch to using the first order Taylor expansion.
    R << 1.0 , -angle_axis(2), angle_axis(1),
        angle_axis(2), 1.0, -angle_axis(0),
        -angle_axis(1), angle_axis(0), 1.0;

    /*
    R[0] =  1.0;
    R[1] = -angle_axis[2];
    R[2] =  angle_axis[1];

    R[3] =  angle_axis[2];
    R[4] =  1.0;
    R[5] = -angle_axis[0];
    R[6] = -angle_axis[1];
    R[7] =  angle_axis[0];
    R[8] = 1.0;*/
  }
}

static Mat3 Hat(const Eigen::Vector3d & v)
{
  Mat3 v_hat;
  v_hat <<
      0.0, -v(2),  v(1),
       v(2), 0.0, -v(0),
      -v(1),  v(0), 0.0;

  return v_hat;
}



// Adjusted from Sophus logAndTheta (so3.hpp)
// https://github.com/strasdat/Sophus/blob/master/sophus/so3.hpp
// Logarithmic map
//
// Computes the logarithm, the inverse of the group exponential which maps
// element of the group (rotation matrices) to elements of the tangent space
// (rotation-vector).
//
// To be specific, this function computes ``vee(logmat(.))`` with
// ``logmat(.)`` being the matrix logarithm and ``vee(.)`` the vee-operator
// of SO(3).
//
static void SO3_logAndTheta(const Mat3 & R, Eigen::Vector3d & omega, double & theta)
{
  using std::sqrt;
  using std::atan;
  using std::abs;
  const double f_epsilon = 1e-6;


  const Eigen::Quaterniond q(R);
  const double squared_n = q.vec().squaredNorm();
  const double n = sqrt(squared_n);
  const double w = q.w();

  double two_atan_nbyw_by_n;

  // Atan-based log thanks to
  //
  // C. Hertzberg et al.:
  // "Integrating Generic Sensor Fusion Algorithms with Sound State
  // Representation through Encapsulation of Manifolds"
  // Information Fusion, 2011
  if (n < f_epsilon) {
    // If quaternion is normalized and n=0, then w should be 1;
    // w=0 should never happen here!
    assert(abs(w) >= f_epsilon);
    double squared_w = w * w;
    two_atan_nbyw_by_n =
        2.0f / w - 2.0f * (squared_n) / (w * squared_w);
  } else {
    if (abs(w) < f_epsilon) {
      if (w > 0.0f) {
        two_atan_nbyw_by_n = M_PI / n;
      } else {
        two_atan_nbyw_by_n = -M_PI / n;
      }
    } else {
      two_atan_nbyw_by_n = 2.0f * atan(n / w) / n;
    }
  }

  theta = two_atan_nbyw_by_n * n;
  omega = two_atan_nbyw_by_n* q.vec();
}

// Taken from Viorela Ila's matlab code (IntrinsicSLAM)
static void SO3_expAndTheta(const Eigen::Vector3d & omega, Mat3 & R, double & theta)
{
  const double f_epsilon = 1e-6;

  theta = omega.norm();
  const double sin_theta = sin(theta);
  const double cos_theta = cos(theta);
  const  Mat3 I = Mat3::Identity();

  const Mat3 Omega = Hat(omega);
  const Mat3 Omega2 = Omega * Omega;

  if (abs(theta) < f_epsilon)
  {
    R = I + Omega + 0.5 * Omega2;
  }
  else
  {
    R = I + (sin_theta/theta) * Omega + ((1-cos_theta)/(theta*theta)) * Omega2;
  }
}


// Adjusted from Sophus log (sim3.hpp line 440)
// https://github.com/strasdat/Sophus/blob/master/sophus/sim3.hpp
// Logarithmic map
//
// Computes the logarithm, the inverse of the group exponential which maps
// element of the group (rigid body transformations) to elements of the
// tangent space (twist).
//
// To be specific, this function computes ``vee(logmat(.))`` with
// ``logmat(.)`` being the matrix logarithm and ``vee(.)`` the vee-operator
// of Sim(3).
//
static void Sim3_log(const Mat3 & R, const Vec3 & t, const double & scale, Eigen::Matrix<double, 7, 1> & v_log)
{
  // The derivation of the closed-form Sim(3) logarithm for is done
  // analogously to the closed-form solution of the SE(3) logarithm, see
  // J. Gallier, D. Xu, "Computing exponentials of skew symmetric matrices and
  // logarithms of orthogonal matrices", IJRA 2002.
  // https://pdfs.semanticscholar.org/cfe3/e4b39de63c8cabd89bf3feff7f5449fc981d.pdf
  // (Sec. 6., pp. 8)
  const double f_epsilon = 1e-6;

  // Sigma
  double sigma = log(scale);

  // Theta and omega
  double theta;
  Vec3 omega;
  SO3_logAndTheta(R,omega,theta);

  // upsilon
  using std::abs;
  using std::sin;
  using std::abs;
  const  Mat3 I = Mat3::Identity();

  const Mat3 Omega = Hat(omega);
  const Mat3 Omega2 = Omega * Omega;
  const double scale_sq = scale * scale;
  const double theta_sq = theta * theta;
  const double sin_theta = sin(theta);
  const double cos_theta = cos(theta);

  double a, b, c;
  if (abs(sigma * sigma) < f_epsilon) {
    c = 1.0f - (0.5f * sigma);
    a = -0.5f;
    if (abs(theta_sq) <f_epsilon) {
      b = 1.0f / 12.0f;
    } else {
      b = (theta * sin_theta + 2.0f * cos_theta - 2.0f) /
          (2.0f * theta_sq * (cos_theta - 1.0f));
    }
  } else {
    const double scale_cu = scale_sq * scale;
    c = sigma / (scale - 1.0f);
    if (abs(theta_sq) < f_epsilon) {
      a = (-sigma * scale + scale - 1.0f) / ((scale - 1.0f) * (scale - 1.0f));
      b = (scale_sq * sigma - 2.0f * scale_sq + scale * sigma + 2.0f * scale) /
          (2.0f * scale_cu - 6.0f * scale_sq + 6.0f * scale - 2.0f);
    } else {
      const double s_sin_theta = scale * sin_theta;
      const double s_cos_theta = scale * cos_theta;
      a = (theta * s_cos_theta - theta - sigma * s_sin_theta) /
          (theta * (scale_sq - 2.0f * s_cos_theta + 1.0f));
      b = -scale *
          (theta * s_sin_theta - theta * sin_theta + sigma * s_cos_theta -
           scale * sigma + sigma * cos_theta - sigma) /
          (theta_sq * (scale_cu - 2.0f * scale * s_cos_theta - scale_sq +
                       2.0f * s_cos_theta + scale - 1.0f));
    }
  }
  const Mat3 W_inv = a * Omega + b * Omega2 + c * I;
  // upsilon
  v_log.template head<3>() = W_inv * t;
  // omega
  v_log.template tail<4>().template head<3>() = omega;
  // sigma
  v_log(6) = sigma;

}


// Adjusted from Sophus exp (sim3.hpp line 280)
// https://github.com/strasdat/Sophus/blob/master/sophus/sim3.hpp
// Group exponential
//
// This functions takes in an element of tangent space and returns the
// corresponding element of the group Sim(3).
//
// The first three components of ``a`` represent the translational part
// ``upsilon`` in the tangent space of Sim(3), the following three components
// of ``a`` represents the rotation vector ``omega`` and the final component
// represents the logarithm of the scaling factor ``sigma``.
// To be more specific, this function computes ``expmat(hat(a))`` with
// ``expmat(.)`` being the matrix exponential and ``hat(.)`` the hat-operator
// of Sim(3), see below.
//
static void Sim3_exp(const Eigen::Matrix<double, 7, 1> & v_log, Mat3 & R, Vec3 & t, double & scale)
{
  const double f_epsilon = 1e-6;
  // For the derivation of the exponential map of Sim(3) see
  // H. Strasdat, "Local Accuracy and Global Consistency for Efficient Visual
  // SLAM", PhD thesis, 2012.
  // http://hauke.strasdat.net/files/strasdat_thesis_2012.pdf (A.5, pp. 186)
  const Eigen::Vector3d upsilon = v_log.template head<3>();
  const Eigen::Vector3d omega = v_log.template tail<4>().template head<3>();
  const double sigma = v_log(6);
  // Scale
  scale = exp(sigma);

  const double theta = omega.norm();
  const Mat3 Omega = Hat(omega);
  const Mat3 Omega2 = Omega * Omega;


  using std::abs;
  using std::sin;
  using std::cos;
  const Mat3 I = Mat3::Identity();
  double A, B, C;

  if (abs(sigma) < f_epsilon) {
    C = 1.0;
    if (abs(theta) < f_epsilon) {
      A = 0.5;
      B = 1.0 / 6.0;
    } else {
      const double theta_sq = theta * theta;
      A = (1.0 - cos(theta)) / theta_sq;
      B = (theta - sin(theta)) / (theta_sq * theta);
    }
  } else {
    C = (scale - 1.0) / sigma;
    if (abs(theta) < f_epsilon) {
      const double sigma_sq = sigma * sigma;
      A = ((sigma - 1.0) * scale + 1.0) / sigma_sq;
      B = ((0.5 * sigma * sigma - sigma + 1.0) * scale) / (sigma_sq * sigma);
    } else {
      const double theta_sq = theta * theta;
      const double a = scale * sin(theta);
      const double b = scale * cos(theta);
      const double c = theta_sq + sigma * sigma;
      A = (a * sigma + (1.0 - b) * theta) / (theta * c);
      B = (C - ((b - 1.0) * sigma + a * theta) / (c)) * 1.0 / (theta_sq);
    }
  }

  // Translation
  const Mat3 W = A * Omega + B * Omega2 + C * I;
  t = W * upsilon;

  // Rotation
  double th;
  SO3_expAndTheta(omega,R,th);
}

}
}
