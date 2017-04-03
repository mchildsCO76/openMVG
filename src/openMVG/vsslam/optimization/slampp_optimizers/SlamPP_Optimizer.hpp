#pragma once

/*
 * SlamPP_Optimizer - Global cameras and global XYZ landmarks
 */

#define SIM3_USE_ROBUST_EDGES
#define USE_EXPLICIT_HARD_UF // use explicit UF via a const vertex?
#include <string.h>
#include <stdio.h>


#include "slam/Sim3SolverBase.h" // C3DJacobians, CBAJacobians, generally useful functions for BA and SE(3), does not need to be included
#define LANDMARK_TYPE_XYZ
#define LANDMARKS_GLOBAL

#define SIM3_USE_ROBUST_EDGES
#define USE_EXPLICIT_HARD_UF // use explicit UF via a const vertex?

#include "slam/LinearSolver_UberBlock.h"
//#include "slam/LinearSolver_CholMod.h" // not used
#include "slam/ConfigSolvers.h" // nonlinear graph solvers
//#include "slam/Timer.h" // not used
#include "slam/3DSolverBase.h" // want C3DJacobians::Quat_to_AxisAngle() and C3DJacobians::AxisAngle_to_Quat()
#include "slam/Sim3_Types.h"
//#include "slam/Marginals.h" // included from nonlinear solver, if supported
#include "slam/NonlinearSolver_Lambda_DL.h"
//#include "slam/BASolverBase.h" // included from BA_Types.h
#include "slam/Eigenvalues.h"



#include <map>


/**
 *	@brief bundle adjustment optimizer
 *	@note The purpose is not to have any interaction with CSystemType or CNonlinearSolverType
 *		in the header, in order to move the slow building of SLAM++ to BAOptimizer.cpp, which
 *		is rarely changed throughout the development process.
 */
namespace openMVG {
namespace VSSLAM {

template <class CSolverType, const bool b_can_delay>
class CDelayedOptimizationCaller {
public:
  static void Delay(CSolverType &UNUSED(r_solver))
  {}

  static void Enable(CSolverType &UNUSED(r_solver))
  {}
};

template <class CSolverType>
class CDelayedOptimizationCaller<CSolverType, true> {
public:
  static void Delay(CSolverType &r_solver)
  {
    r_solver.Delay_Optimization();
  }

  static void Enable(CSolverType &r_solver)
  {
    r_solver.Enable_Optimization();
  }
};

class SlamPP_Optimizer {
public:
  size_t m_undefined_camera_id; // index of undefined camera (if we dont have constant vertices its size_t(-1) otherwise size_t::max/2 +1

  SlamPP_Optimizer(size_t undefined_camera_id) : m_undefined_camera_id(undefined_camera_id){};

	virtual ~SlamPP_Optimizer() = default;

	virtual size_t n_Vertex_Num() const = 0;


	virtual Eigen::Map<const Eigen::VectorXd> r_Vertex_State(size_t n_index) const = 0; // returns const map of vertex state only
	virtual Eigen::Map<Eigen::VectorXd> r_Vertex_State(size_t n_index) = 0; // returns map of vertex state only


  virtual Eigen::Vector3d v_XYZVertex(size_t n_index) const = 0;

	virtual void Delay_Optimization() = 0;
	virtual void Enable_Optimization() = 0; // throw(srd::bad_alloc, std::runtime_error)

	virtual void Set_TrustRadius(double f_trust_radius) = 0;
	virtual void Set_TrustRadius_Persistence(bool b_trust_radius_persistent) = 0;
	virtual void Set_UpdateThreshold(double f_update_thresh) = 0;
	virtual void Set_AllBatch(bool b_all_batch) = 0;
	virtual void Optimize(size_t n_max_iteration_num = 5, double f_min_dx_norm = .01, double f_min_dl_norm = .01) = 0; // throw(srd::bad_alloc, std::runtime_error)

	virtual double * Add_CamVertex(size_t n_vertex_id, const Eigen::Matrix<double, 12, 1> &v_cam_state) = 0; // throw(srd::bad_alloc)
  virtual double * Add_CamVertexConst(size_t n_vertex_id, const Eigen::Matrix<double, 12, 1> &v_cam_state) = 0; // throw(srd::bad_alloc)
  virtual double * Add_CamVertexFixed(size_t n_vertex_id, const Eigen::Matrix<double, 12, 1> &v_cam_state) =0; // throw(srd::bad_alloc)

	// a vertex can be added either way, it gets automatically converted to the representation space
	virtual double * Add_XYZVertex(size_t n_vertex_id, size_t n_owner_id, const Eigen::Vector3d &v_xyz_position) = 0; // throw(srd::bad_alloc)
	virtual double * Add_XYZVertexFixed(size_t n_vertex_id, size_t n_owner_id, const Eigen::Vector3d &v_xyz_position) = 0; // throw(srd::bad_alloc)



	virtual void Add_P2CSim3GEdge(size_t n_landmark_vertex_id, size_t n_cam_vertex_id,
		const Eigen::Vector2d &v_observation, const Eigen::Matrix2d &t_information) = 0; // throw(srd::bad_alloc)

	virtual void Show_Stats(bool b_calculate_eigenvalues = true) const = 0;

	virtual bool Dump_State(const char *p_s_filename) const = 0; // note that this saves the state as Sim(3)
	virtual bool Dump_Marginals(const char *p_s_filename) const = 0; // note that this saves the marginals of (relative?) Sim(3)

	virtual bool Dump_Poses_SE3(const char *p_s_filename) const = 0; // only poses, not landmarks
	virtual bool Dump_State_SE3(const char *p_s_filename) const = 0; // note that this converts the state to SE(3) for viewing by our utility which does not support Sim(3) now
	virtual bool Dump_Graph_SE3(const char *p_s_filename) const = 0; // note that this converts the graph to SE(3) for viewing by our utility which does not support Sim(3) now
};

}
}
