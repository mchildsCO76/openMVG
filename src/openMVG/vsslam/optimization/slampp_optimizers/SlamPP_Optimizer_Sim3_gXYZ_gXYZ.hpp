#pragma once

#define SIM3_USE_ROBUST_EDGES
#define USE_EXPLICIT_HARD_UF // use explicit UF via a const vertex?

/*
 * SlamPP_Optimizer_Sim3_gXYZ_gXYZ - Global cameras and global XYZ landmarks
 */

//#include <openMVG/vsslam/optimization/slampp_optimizers/CBAOptimizerCore_Sim3_gXYZ_gXYZ.app>
#include <openMVG/vsslam/optimization/slampp_optimizers/SlamPP_Optimizer.hpp>
#include <map>

/*#include "slam/LinearSolver_UberBlock.h"
#include "slam/ConfigSolvers.h" // nonlinear graph solvers
#include "slam/Sim3_Types.h"*/

#include "slam/Sim3SolverBase.h" // C3DJacobians, CBAJacobians, generally useful functions for BA and SE(3), does not need to be included
#define LANDMARK_TYPE_XYZ
#define LANDMARKS_GLOBAL

#define SIM3_USE_ROBUST_EDGES
#define USE_EXPLICIT_HARD_UF // use explicit UF via a const vertex?
#include <map>


/**
 *	@brief bundle adjustment optimizer
 *	@note The purpose is not to have any interaction with CSystemType or CNonlinearSolverType
 *		in the header, in order to move the slow building of SLAM++ to BAOptimizer.cpp, which
 *		is rarely changed throughout the development process.
 */
namespace openMVG {
namespace VSSLAM {


class SlamPP_Optimizer_Sim3_gXYZ_gXYZ : public SlamPP_Optimizer {
public:
  //int n_dummy_param=1;
  class CBAOptimizerCore_Sim3_gXYZ_gXYZ;
protected:

  CBAOptimizerCore_Sim3_gXYZ_gXYZ * m_p_optimizer; // PIMPL
	std::map<size_t, size_t> m_camera_ownerships; // using a side channel now. this should be a part of the CVertexInvDepth class but I actually want to experiment a bit with global / local convergence so I'm leaving it here for now

public:
	SlamPP_Optimizer_Sim3_gXYZ_gXYZ(size_t undefined_camera_id = size_t(-1), bool b_verbose = false, bool b_use_schur = true,
		bool b_do_marginals = false, bool b_do_icra_style_marginals = false); // throw(srd::bad_alloc) // todo - add more optimization settings as required
	~SlamPP_Optimizer_Sim3_gXYZ_gXYZ();

	size_t n_Vertex_Num() const override;
/*
	CBAOptimizerCore::CSystemType::_TyConstVertexRef r_Vertex(size_t n_index) const override;
	CBAOptimizerCore::CSystemType::_TyVertexRef r_Vertex(size_t n_index) override;*/

	Eigen::Map<const Eigen::VectorXd> r_Vertex_State(size_t n_index) const override; // returns const map of vertex state only
	Eigen::Map<Eigen::VectorXd> r_Vertex_State(size_t n_index) override; // returns map of vertex state only

  Eigen::Vector3d v_XYZVertex(size_t n_index) const override;

	void Delay_Optimization() override;
	void Enable_Optimization() override; // throw(srd::bad_alloc, std::runtime_error)

	void Set_TrustRadius(double f_trust_radius) override;
	void Set_TrustRadius_Persistence(bool b_trust_radius_persistent) override;
	void Set_UpdateThreshold(double f_update_thresh) override;
	void Set_AllBatch(bool b_all_batch) override;
	void Optimize(size_t n_max_iteration_num = 5, double f_min_dx_norm = .01, double f_min_dl_norm = .01) override; // throw(srd::bad_alloc, std::runtime_error)

	double * Add_CamVertex(size_t n_vertex_id, const Eigen::Matrix<double, 12, 1> &v_cam_state) override; // throw(srd::bad_alloc)
	double * Add_CamVertexFixed(size_t n_vertex_id, const Eigen::Matrix<double, 12, 1> &v_cam_state) override; // throw(srd::bad_alloc)
  double * Add_CamVertexConst(size_t n_vertex_id, const Eigen::Matrix<double, 12, 1> &v_cam_state) override; // throw(srd::bad_alloc)

  // a vertex can be added either way, it gets automatically converted to the representation space
  double * Add_XYZVertex(size_t n_vertex_id, size_t n_owner_id, const Eigen::Vector3d &v_xyz_position) override;
  double * Add_XYZVertexFixed(size_t n_vertex_id, size_t n_owner_id, const Eigen::Vector3d &v_xyz_position) override;

  // Observations
	void Add_P2CSim3GEdge(size_t n_landmark_vertex_id, size_t n_cam_vertex_id,
		const Eigen::Vector2d &v_observation, const Eigen::Matrix2d &t_information) override; // throw(srd::bad_alloc)

	void Show_Stats(bool b_calculate_eigenvalues = true) const override;

	bool Dump_State(const char *p_s_filename) const override; // note that this saves the state as Sim(3)
	bool Dump_Marginals(const char *p_s_filename) const override; // note that this saves the marginals of (relative?) Sim(3)

	bool Dump_Poses_SE3(const char *p_s_filename) const override; // only poses, not landmarks
	bool Dump_State_SE3(const char *p_s_filename) const override; // note that this converts the state to SE(3) for viewing by our utility which does not support Sim(3) now
	bool Dump_Graph_SE3(const char *p_s_filename) const override; // note that this converts the graph to SE(3) for viewing by our utility which does not support Sim(3) now


protected:
	Eigen::Vector3d v_LandmarkState_XYZ(size_t n_index) const;

private:
	SlamPP_Optimizer_Sim3_gXYZ_gXYZ(const SlamPP_Optimizer_Sim3_gXYZ_gXYZ &r_other); // no copy
	SlamPP_Optimizer_Sim3_gXYZ_gXYZ &operator =(const SlamPP_Optimizer_Sim3_gXYZ_gXYZ &r_other); // no copy
};

}
}
