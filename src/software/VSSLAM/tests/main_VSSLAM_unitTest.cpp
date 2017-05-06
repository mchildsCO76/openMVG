// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Klemen Istenic

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.




#include <openMVG/types.hpp>
#include <openMVG/vsslam/optimization/Sim3.hpp>

#include "slam/3DSolverBase.h" // want C3DJacobians::Quat_to_AxisAngle() and C3DJacobians::AxisAngle_to_Quat()
#include "slam/Sim3SolverBase.h"

using namespace openMVG;
using namespace openMVG::vsslam;

bool check_equal (Mat4 &T_A, Mat4 &T_B)
{
  bool b_equal = true;
  for (size_t i=0; i < 4; i++)
  {
    for (size_t j=0; j < 4; j++)
    {
      if (fabs(T_A(i,j) - T_B(i,j)) > 0.001)
      {
        std::cout<<"Error: T_1: "<<T_A(i,j)<<" T_2: "<<T_B(i,j)<<"\n";
        b_equal = false;
      }
    }
  }

  if (!b_equal)
  {
    std::cout<<"Matrices are not equal: \nT_1:\n"<<T_A<<"\nT_2:\n"<<T_B<<"\n";
  }
  return b_equal;
}

bool check_equal (Vec7 &v_A, Vec7 &v_B)
{
  bool b_equal = true;
  for (size_t i=0; i < 7; i++)
  {
    if (fabs(v_A(i) - v_B(i)) > 0.001)
    {
      std::cout<<"Error: v_1: "<<v_A(i)<<" v_2: "<<v_B(i)<<"\n";
      b_equal = false;
    }
  }

  if (!b_equal)
  {
    std::cout<<"Vectors are not equal: \nv_1:\n"<<v_A<<"v_2:\n"<<v_B<<"\n";
  }
  return b_equal;
}

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << "VISUAL S-SLAM -- Unit test --" << std::endl;

  Mat4 T_A_gt, T_B_gt, T_C_gt, T_D_gt;
  Mat4 T_A, T_B, T_C, T_D;
  Vec7 v_A, v_B, v_C, v_D;  // v,omega,sigma

  /////////////////////////////
  // Dataset - translation
  /////////////////////////////

  std::cout<<" Test A:\n";
  double s_A = 1.0;
  Mat3 R_A = Mat3::Identity();
  Vec3 t_A = Vec3(5.0,1.0,5.0);
  /*T_A_gt = Mat4::Identity();
  T_A_gt.block(0,0,3,3) = s_A*R_A;
  T_A_gt.block(0,3,3,1) = t_A;*/
  T_A_gt << 0.6928, 0.0, 0.4, 2.0,  0.0, 0.8, 0.0, 2.0, -0.4, 0.0, 0.6928, 2.0, 0.0,0.0,0.0, 1.0;
  Vec7 v_A_gt;
  //v_A_gt << 5.0,1.0,5.0,0.0,0.0,0.0,0.0;
  v_A_gt << 2.4777, 2.0, 1.4305, 0.0, 0.5236, 0.0, -0.2231;

  //std::cout<<"AAA: "<<Sophus::Sim3d::exp(v_A_gt).matrix()<<"\n";
  //std::cout<<"BBB: "<<Sophus::Sim3d(T_A_gt).log()<<"\n";
  vsslam::Sim3_log(T_A_gt, v_A);
  vsslam::Sim3_exp(v_A_gt, T_A);

  bool b_test_A;
  b_test_A = check_equal(T_A_gt,T_A);
  b_test_A = b_test_A && check_equal(v_A_gt,v_A);

  std::cout<<" Test A:\n"
      <<"T_A gt:\n"<<T_A_gt<<"\n"
      <<"T_A:\n"<<T_A<<"\n"
      <<"v_A gt:\n"<<v_A_gt.transpose()<<"\n"
      <<"v_A:\n"<<v_A.transpose()<<"\n";
  if (b_test_A)
    std::cout<<"Test A: SUCCESSFUL\n";
  else
    std::cout<<"Test A: FAILED\n";

  /*
  // SLAM++
  Eigen::Quaternionf q_rot(T_A_gt.block(0,0,3,3));
  double q_scale = q_rot.norm();
  q_rot.normalize();
  double T_scale = T_A_gt.block(0,0,3,1).norm();
  Vec3 T_translation = T_A_gt.block(0,3,3,1);

  std::cout<<"Scale A: "<<q_scale<<" :: "<<T_scale<<"\n";

  CSim3Jacobians::TSim3 pose (T_translation, q_rot, T_scale);// = CSim3Jacobians::t_Exp(vertex); // not Sim3 yet
*/
  Eigen::Matrix4d T;
  T << 0.6928, 0, 0.4000, 2.0000,
      0, 0.8000, 0, 2.0000,
      -0.4000, 0, 0.6928, 2.0000,
      0, 0, 0, 1.0000;

  Eigen::Vector7d v_tRs;
  v_tRs.head<3>() = T.topRightCorner<3, 1>(); // translation
  v_tRs(6) = T.topLeftCorner<3, 3>().norm() / sqrt(3.0); // scale

  Eigen::Vector3d v_temp;
  C3DJacobians::Quat_to_AxisAngle(Eigen::Quaterniond(T.topLeftCorner<3, 3>() / v_tRs(6)).normalized(),
         v_temp); // rotation
   v_tRs.segment<3>(3) = v_temp;


  std::cout << "T =" << std::endl << T << std::endl << std::endl;
  std::cout << "v_tRs =" << std::endl << v_tRs << std::endl << std::endl;

  CSim3Jacobians::TSim3 t_sim3(v_tRs, CSim3Jacobians::TSim3::from_tRs_vector);

  Eigen::Vector7d v_log = t_sim3.v_Log();

  std::cout << "v_log =" << std::endl << v_log << std::endl;

  CSim3Jacobians::TSim3 t_sim3B(v_log, CSim3Jacobians::TSim3::from_sim3_vector);

  Mat4 T2 = Mat4::Identity();
  T2.block(0,0,3,3) = t_sim3B.f_scale * t_sim3B.t_rotation.toRotationMatrix();
  T2.block(0,3,3,1) = t_sim3B.v_translation;



  std::cout << "T2 =" << std::endl << T2 << std::endl << std::endl;

  /////////////////////////////
  // Dataset 2 - rotation (eul2rotm([0.8,0.2,0.3]))
  /////////////////////////////

  std::cout<<" Test B:\n";
  double s_B = 1.0;
  Mat3 R_B;
  R_B <<   0.682818960388909, -0.644412239881579, 0.344225409323913, 0.703056729101466, 0.707705892854635, -0.069739550356817, -0.198669330795061,   0.289629477625516,   0.936293363584199;
  Vec3 t_B = Vec3(0.0,0.0,0.0);
  T_B_gt = Mat4::Identity();
  T_B_gt.block(0,0,3,3) = s_B*R_B;
  T_B_gt.block(0,3,3,1) = t_B;
  Vec7 v_B_gt;
  v_B_gt << 0.0,0.0,0.0,0.203019577517389,  0.306699387518686, 0.761230266104565, 0.0;

  vsslam::Sim3_log(T_B_gt, v_B);
  vsslam::Sim3_exp(v_B_gt, T_B);

  bool b_test_B;
  b_test_B = check_equal(T_B_gt,T_B);
  b_test_B = b_test_B && check_equal(v_B_gt,v_B);


  std::cout<<" Test B:\n"
      <<"T_B gt:\n"<<T_B_gt<<"\n"
      <<"T_B:\n"<<T_B<<"\n"
      <<"v_B gt:\n"<<v_B_gt<<"\n"
      <<"v_B:\n"<<v_B<<"\n";
  if (b_test_B)
    std::cout<<"Test B: SUCCESSFUL\n";
  else
    std::cout<<"Test B: FAILED\n";

  /////////////////////////////
  // Dataset 3 rot+trans:  eul2rotm([0.8,0.2,0.3]) t = [-3.2,1.5,-6.4]
  /////////////////////////////
  std::cout<<" Test C:\n";
  double s_C = 1.0;
  Mat3 R_C;
  R_C <<  0.682818960388909, -0.644412239881579, 0.344225409323913, 0.703056729101466, 0.707705892854635, -0.069739550356817, -0.198669330795061,   0.289629477625516,   0.936293363584199;
  Vec3 t_C = Vec3(-3.2, 1.5, -6.4);
  T_C_gt = Mat4::Identity();
  T_C_gt.block(0,0,3,3) = s_C*R_C;
  T_C_gt.block(0,3,3,1) = t_C;
  Vec7 v_C_gt;
  v_C_gt <<  -4.646119579001327, 0.710336676023224, -5.696165867760540,  0.203019577517389,  0.306699387518686, 0.761230266104565, 0.0;


  vsslam::Sim3_log(T_C_gt, v_C);
  vsslam::Sim3_exp(v_C_gt, T_C);

  bool b_test_C;
  b_test_C = check_equal(T_C_gt,T_C);
  b_test_C = b_test_C && check_equal(v_C_gt,v_C);

  std::cout<<" Test C:\n"
      <<"T_C gt:\n"<<T_C_gt<<"\n"
      <<"T_C:\n"<<T_C<<"\n"
      <<"v_C gt:\n"<<v_C_gt<<"\n"
      <<"v_C:\n"<<v_C<<"\n";
  if (b_test_C)
    std::cout<<"Test C: SUCCESSFUL\n";
  else
    std::cout<<"Test C: FAILED\n";


  /////////////////////////////
  // Dataset 4 rot+trans+scale:  eul2rotm([0.8,0.2,0.3]) t = [-3.2,1.5,-6.4]
  /////////////////////////////
  std::cout<<" Test D:\n";
  double s_D = 1.4;
  Mat3 R_D;
  R_D << 0.682818960388909, -0.644412239881579, 0.344225409323913, 0.703056729101466, 0.707705892854635, -0.069739550356817, -0.198669330795061,   0.289629477625516,   0.936293363584199;
  Vec3 t_D = Vec3(-3.2, 1.5, -6.4);
  T_D_gt = Mat4::Identity();
  T_D_gt.block(0,0,3,3) = s_D*R_D;
  T_D_gt.block(0,3,3,1) = t_D;
  Vec7 v_D_gt;
  v_D_gt << -4.646119579001327, 0.710336676023224, -5.696165867760540,  0.203019577517389,  0.306699387518686, 0.761230266104565, 0.3365;


  //T_D = Sophus::Sim3d::exp(v_D_gt).matrix();
  //v_D = Sophus::Sim3d(T_D_gt).log();
  vsslam::Sim3_log(T_D_gt, v_D);
  vsslam::Sim3_exp(v_D_gt, T_D);

  Mat4 T_DD;
  vsslam::Sim3_exp(v_D, T_DD);

  bool b_test_D;
  b_test_D = check_equal(T_D_gt,T_D);
  b_test_D = b_test_D && check_equal(v_D_gt,v_D);

  std::cout<<" Test D:\n"
      <<"T_D gt:\n"<<T_D_gt<<"\n"
      <<"T_D:\n"<<T_D<<"\n"
      <<"T_DD:\n"<<T_DD<<"\n"
      <<"v_D gt:\n"<<v_D_gt<<"\n"
      <<"v_D:\n"<<v_D<<"\n";

  if (b_test_D)
    std::cout<<"Test D: SUCCESSFUL\n";
  else
    std::cout<<"Test D: FAILED\n";


}
