#include <openMVG/vsslam/optimization/slampp_optimizers/SlamPP_Optimizer_Sim3_gXYZ_gXYZ.hpp>

#include <openMVG/vsslam/optimization/slampp_optimizers/SlamPP_IO.hpp>
#include <openMVG/vsslam/optimization/slampp_optimizers/SlamPP_Utils.hpp>

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

namespace openMVG  {
namespace VSSLAM  {



class SlamPP_Optimizer_Sim3_gXYZ_gXYZ::CBAOptimizerCore_Sim3_gXYZ_gXYZ {
public:
  // landmark: g_xyz
  using _TyLandmark = CVertexXYZ;
  using _TyObservation = CEdgeP2C_XYZ_Sim3_G;

  /*
  // landmark l_xyz
  using _TyLandmark = CVertexXYZ;
  using _TyObservationS = CEdgeP2C_XYZ_Sim3_LS;
  using _TyObservationO = CEdgeP2C_XYZ_Sim3_LO;

  // landmark l_invdepth
  using _TyLandmark = CVertexInvDepth;
  using _TyObservationS = CEdgeP2C_XYZ_Sim3_LS;
  using _TyObservationO = CEdgeP2C_XYZ_Sim3_LO;
*/
  using TVertexTypelist = MakeTypelist_Safe((CVertexCamSim3, _TyLandmark));
  using TConstVertexTypelist = MakeTypelist_Safe((CVertexCamSim3, _TyLandmark));
  using TEdgeTypelist = MakeTypelist_Safe((_TyObservation, CEdgePoseCamSim3));
  using CSystemType = CFlatSystem<CBaseVertex, TVertexTypelist, CBaseEdge, TEdgeTypelist, CNullUnaryFactorFactory,
      CBaseVertex, TConstVertexTypelist>;

  using CLinearSolverType = CLinearSolver_UberBlock<CSystemType::_TyHessianMatrixBlockList>;
  using CNonlinearSolverType = CNonlinearSolver_Lambda_DL<CSystemType, CLinearSolverType>;

protected:
  CSystemType m_system;
  CNonlinearSolverType m_solver;

public:
  inline CBAOptimizerCore_Sim3_gXYZ_gXYZ(bool b_verbose = false, bool b_use_schur = true,
    bool b_do_marginals = false, bool b_do_icra_style_marginals = false)
    :m_solver(m_system, TIncrementalSolveSetting(),
    TMarginalsComputationPolicy((b_do_marginals)? marginals::do_calculate : marginals::do_not_calculate,
    frequency::Never(), mpart_Diagonal, mpart_Diagonal), // batch marginals
    b_verbose, CLinearSolverType(), b_use_schur)
  {
    m_solver.Set_ICRA15_Style_IncrementalMargs(b_do_icra_style_marginals); // !!
  }

  CSystemType::_TyConstVertexRef r_Vertex(size_t n_index) const
  {
    return m_system.r_Vertex_Pool()[n_index];
  }

  CSystemType::_TyVertexRef r_Vertex(size_t n_index)
  {
    return m_system.r_Vertex_Pool()[n_index];
  }

  inline CSystemType &r_System()
  {
    return m_system;
  }

  inline CNonlinearSolverType &r_Solver()
  {
    return m_solver;
  }

  inline const CSystemType &r_System() const
  {
    return m_system;
  }

  inline const CNonlinearSolverType &r_Solver() const
  {
    return m_solver;
  }

private:
  CBAOptimizerCore_Sim3_gXYZ_gXYZ(const CBAOptimizerCore_Sim3_gXYZ_gXYZ &r_other); // no copy
  CBAOptimizerCore_Sim3_gXYZ_gXYZ &operator =(const CBAOptimizerCore_Sim3_gXYZ_gXYZ &r_other); // no copy
};


/*
 *								=== SlamPP_Optimizer_Sim3_gXYZ_gXYZ_Sim3_gXYZ_gXYZ ===
 */

SlamPP_Optimizer_Sim3_gXYZ_gXYZ::SlamPP_Optimizer_Sim3_gXYZ_gXYZ
(
  size_t undefined_camera_id,
  bool b_verbose,
  bool b_use_schur,
	bool b_do_marginals,
	bool b_do_icra_style_marginals) // throw(srd::bad_alloc)
	:SlamPP_Optimizer(undefined_camera_id), m_p_optimizer(new CBAOptimizerCore_Sim3_gXYZ_gXYZ(b_verbose, b_use_schur,
	 b_do_marginals, b_do_icra_style_marginals))
{}

SlamPP_Optimizer_Sim3_gXYZ_gXYZ::~SlamPP_Optimizer_Sim3_gXYZ_gXYZ()
{
	delete m_p_optimizer;
}

size_t SlamPP_Optimizer_Sim3_gXYZ_gXYZ::n_Vertex_Num() const
{
	return m_p_optimizer->r_System().r_Vertex_Pool().n_Size();
}


Eigen::Vector3d SlamPP_Optimizer_Sim3_gXYZ_gXYZ::v_LandmarkState_XYZ(size_t n_vertex_id) const
{
  return m_p_optimizer->r_System().r_Vertex_Pool().r_At<CBAOptimizerCore_Sim3_gXYZ_gXYZ::_TyLandmark>(n_vertex_id).r_v_State();
}

Eigen::Map<const Eigen::VectorXd> SlamPP_Optimizer_Sim3_gXYZ_gXYZ::r_Vertex_State(size_t n_index) const
{
	return ((CBAOptimizerCore_Sim3_gXYZ_gXYZ::CSystemType::_TyConstVertexRef)m_p_optimizer->r_Vertex(n_index)).v_StateC(); // explicit const required in windows for some reason // todo - this was working before; investigate
}

Eigen::Map<Eigen::VectorXd> SlamPP_Optimizer_Sim3_gXYZ_gXYZ::r_Vertex_State(size_t n_index)
{
	return m_p_optimizer->r_Vertex(n_index).v_State();
}


Eigen::Vector3d SlamPP_Optimizer_Sim3_gXYZ_gXYZ::v_XYZVertex(size_t n_index) const
{
  return v_LandmarkState_XYZ(n_index); // now the vertices are global

}
void SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Delay_Optimization()
{
	CDelayedOptimizationCaller<CBAOptimizerCore_Sim3_gXYZ_gXYZ::CNonlinearSolverType,
	SlamPP_Optimizer_Sim3_gXYZ_gXYZ::CBAOptimizerCore_Sim3_gXYZ_gXYZ::CNonlinearSolverType::solver_HasDelayedOptimization>::Delay(m_p_optimizer->r_Solver());
	// some solvers do not implement this
}

void SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Enable_Optimization() // throw(srd::bad_alloc, std::runtime_error)
{
	CDelayedOptimizationCaller<CBAOptimizerCore_Sim3_gXYZ_gXYZ::CNonlinearSolverType,
	SlamPP_Optimizer_Sim3_gXYZ_gXYZ::CBAOptimizerCore_Sim3_gXYZ_gXYZ::CNonlinearSolverType::solver_HasDelayedOptimization>::Enable(m_p_optimizer->r_Solver());
	// some solvers do not implement this
}

void SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Set_TrustRadius(double f_trust_radius)
{
	m_p_optimizer->r_Solver().Set_StepSize(f_trust_radius);
}

void SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Set_TrustRadius_Persistence(bool b_trust_radius_persistent)
{
	m_p_optimizer->r_Solver().Set_StepSize_Persistence(b_trust_radius_persistent);
}

void SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Set_UpdateThreshold(double f_update_thresh)
{
	m_p_optimizer->r_Solver().Set_UpdateThreshold(f_update_thresh);
}

void SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Set_AllBatch(bool b_all_batch)
{
	m_p_optimizer->r_Solver().Set_AllBatch(b_all_batch);
}

void SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Optimize(size_t n_max_iteration_num /*= 5*/,
	double f_min_dx_norm /*= .01*/, double f_min_dl_norm /*= .01*/) // throw(srd::bad_alloc, std::runtime_error)
{
	m_p_optimizer->r_Solver().Optimize(n_max_iteration_num, f_min_dx_norm, f_min_dl_norm);
}


double * SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Add_CamVertex(size_t n_vertex_id, const Eigen::Matrix<double, 12, 1> &v_cam_state) // throw(srd::bad_alloc)
{
	CVertexCamSim3 &r_cam0 = m_p_optimizer->r_System().r_Get_Vertex<CVertexCamSim3>(n_vertex_id, v_cam_state);
	return r_cam0.r_v_State().data();
}

// Fixed vertex is vertex with unary factor
double * SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Add_CamVertexFixed(size_t n_vertex_id, const Eigen::Matrix<double, 12, 1> &v_cam_state) // throw(srd::bad_alloc)
{
  std::cout<<"AA1\n";
  CVertexCamSim3 &r_cam0 = m_p_optimizer->r_System().r_Get_Vertex<CVertexCamSim3>(n_vertex_id, v_cam_state);
  std::cout<<"AA2\n";
  // Add unary factor (UF) to make the optimized system positive definite
  m_p_optimizer->r_System().r_Add_Edge(CEdgePoseCamSim3(n_vertex_id, size_t(-1),
        Eigen::Vector7d::Zero(), Eigen::Matrix7d::Identity() * 100, m_p_optimizer->r_System()));
  std::cout<<"AA3\n";
  // Return pointer to state of vertex
  return r_cam0.r_v_State().data();
}

// Const vertex has id starting with size_t_max - 1
double * SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Add_CamVertexConst(size_t n_vertex_id, const Eigen::Matrix<double, 12, 1> &v_cam_state) // throw(srd::bad_alloc)
{
  CVertexCamSim3 &r_cam0 = m_p_optimizer->r_System().r_Get_ConstVertex<CVertexCamSim3>(n_vertex_id, v_cam_state);
  // Add unary factor (UF) to make the optimized system positive definite
  m_p_optimizer->r_System().r_Add_Edge(CEdgePoseCamSim3(n_vertex_id, m_undefined_camera_id,
        Eigen::Vector7d::Zero(), Eigen::Matrix7d::Identity() * 100, m_p_optimizer->r_System()));
  // Return pointer to state of vertex
  return r_cam0.r_v_State().data();
}

double * SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Add_XYZVertex(size_t n_vertex_id, size_t n_owner_id, const Eigen::Vector3d &v_xyz_position)
{
  CBAOptimizerCore_Sim3_gXYZ_gXYZ::_TyLandmark &r_landmark = m_p_optimizer->r_System().r_Get_Vertex<CBAOptimizerCore_Sim3_gXYZ_gXYZ::_TyLandmark>(n_vertex_id,
    CParserBase::TVertexXYZ(int(n_vertex_id), v_xyz_position(0), v_xyz_position(1), v_xyz_position(2))); // initialize via parse primitive, let the vertex class convert to its representation
  m_camera_ownerships[n_vertex_id] = n_owner_id; // vertex defined in global space
  return r_landmark.r_v_State().data();
}

double * SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Add_XYZVertexFixed(size_t n_vertex_id, size_t n_owner_id, const Eigen::Vector3d &v_xyz_position)
{
  CBAOptimizerCore_Sim3_gXYZ_gXYZ::_TyLandmark &r_landmark = m_p_optimizer->r_System().r_Get_ConstVertex<CBAOptimizerCore_Sim3_gXYZ_gXYZ::_TyLandmark>(n_vertex_id,
    CParserBase::TVertexXYZ(int(n_vertex_id), v_xyz_position(0), v_xyz_position(1), v_xyz_position(2))); // initialize via parse primitive, let the vertex class convert to its representation

  // Save the owner id
  m_camera_ownerships[n_vertex_id] = n_owner_id; // vertex defined in global space

  // Return pointer to state of vertex
  return r_landmark.r_v_State().data();
}


void SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Add_P2CSim3GEdge(size_t n_landmark_vertex_id, size_t n_cam_vertex_id,
	const Eigen::Vector2d &v_observation, const Eigen::Matrix2d &t_information) // throw(srd::bad_alloc)
{
	m_p_optimizer->r_System().r_Add_Edge(CBAOptimizerCore_Sim3_gXYZ_gXYZ::_TyObservation(n_landmark_vertex_id, n_cam_vertex_id,
		v_observation, t_information, m_p_optimizer->r_System()));
}





/**
 *  OUTPUTS
 */
void SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Show_Stats(bool b_calculate_eigenvalues) const
{
	//m_p_optimizer->r_Solver().Dump_ReprojectionErrors(); // debug

	std::vector<int> visible_points(m_p_optimizer->r_System().r_Vertex_Pool().n_Size(), 0); // 0 = camera / visible, 1 = partially visible, 2 = invisible
	double f_AVG_reprojection_err = 0;
	double f_RMSE_reprojection_err = 0;
	double f_AVG_reprojection_err_visible = 0;
	double f_RMSE_reprojection_err_visible = 0;
	const CBAOptimizerCore_Sim3_gXYZ_gXYZ::CSystemType::_TyEdgeMultiPool &edge_pool = m_p_optimizer->r_System().r_Edge_Pool();
	size_t n_visible_obs_num = 0;
#ifdef SIM3_USE_ROBUST_EDGES
	double f_AVG_reprojection_err_robust = 0;
	double f_RMSE_reprojection_err_robust = 0;
	double f_robust_weight_sum = 0;
#endif // SIM3_USE_ROBUST_EDGES
	if(!edge_pool.b_Empty()) {
		for(size_t i = 0, n = edge_pool.n_Size(); i < n; ++ i) {
			bool b_visible;
			double f_weight;
#if 1
			double f_reprojection_err;
			edge_pool.For_Each(i, i + 1, CCalcRepErr(f_reprojection_err, f_weight, b_visible)); // need to derive the edge type
			size_t n_point_id = edge_pool[i].n_Vertex_Id(0); // landmark id // use facade
#else
			double f_reprojection_err = edge_pool.r_At<_TyObservation>(i).f_Reprojection_Error(b_visible); // know the edge type
			Eigen::Vector1d r; r(0) = f_reprojection_err; // can't initialize directly, the constructor is ambiguous
			f_weight = edge_pool.r_At<_TyObservation>(i).f_RobustWeight(r);
			size_t n_point_id = edge_pool.r_At<_TyObservation>(i).n_Vertex_Id(0); // landmark id // know the edge type
#endif
			visible_points[n_point_id] |= (b_visible)? 1 : 2; // need to combine
			f_AVG_reprojection_err += f_reprojection_err;
			f_RMSE_reprojection_err += f_reprojection_err * f_reprojection_err;
#ifdef SIM3_USE_ROBUST_EDGES
			f_AVG_reprojection_err_robust += f_reprojection_err * f_weight;
			f_RMSE_reprojection_err_robust += f_reprojection_err * f_reprojection_err * f_weight;
			f_robust_weight_sum += f_weight;
#endif // SIM3_USE_ROBUST_EDGES
			if(b_visible) {
				++ n_visible_obs_num;
				f_AVG_reprojection_err_visible += f_reprojection_err;
				f_RMSE_reprojection_err_visible += f_reprojection_err * f_reprojection_err;
			}
		}
		f_AVG_reprojection_err /= edge_pool.n_Size();
		f_RMSE_reprojection_err /= edge_pool.n_Size();
		f_RMSE_reprojection_err = sqrt(f_RMSE_reprojection_err);
#ifdef SIM3_USE_ROBUST_EDGES
		f_AVG_reprojection_err_robust /= f_robust_weight_sum;
		f_RMSE_reprojection_err_robust /= f_robust_weight_sum;
		f_RMSE_reprojection_err_robust = sqrt(f_RMSE_reprojection_err_robust);
#endif // SIM3_USE_ROBUST_EDGES
		f_AVG_reprojection_err_visible /= n_visible_obs_num;
		f_RMSE_reprojection_err_visible /= n_visible_obs_num;
		f_RMSE_reprojection_err_visible = sqrt(f_RMSE_reprojection_err_visible);
	}

	size_t n_fully_visible = 0;
	size_t n_fully_invisible = 0;
	size_t n_partly_invisible = 0;
	FILE *p_fw = fopen("visibility_info.txt", "w");
	for(size_t i = 0, n = visible_points.size(); i < n; ++ i) {
		if(visible_points[i] == 0) {
			// a camera, do nothing
		} else if(visible_points[i] == 1) {
			visible_points[i] = 0; // fully visible
			++ n_fully_visible;
		} else if(visible_points[i] == 2)
			++ n_fully_invisible; // fully invisible
		else if(visible_points[i] == 3) {
			visible_points[i] = 1; // visible from some cameras, invisible from the others
			++ n_partly_invisible;
		}

		if(p_fw) {
			double f_value = visible_points[i] * 10;
			size_t n_dim = m_p_optimizer->r_System().r_Vertex_Pool()[i].n_Dimension();
			for(size_t j = 0; j < n_dim; ++ j)
				fprintf(p_fw, (j + 1 < n_dim)? "%f " : "%f\n", f_value);
		}
		// write "fake marginals" that can be displayed by the graph viewer as false colours
	}
	if(p_fw)
		fclose(p_fw);
	// transform the visibility flags

	printf("reprojection error avg: %g px, RMSE: %g px (" PRIsize " observations)\n",
		f_AVG_reprojection_err, f_RMSE_reprojection_err, edge_pool.n_Size());
#ifdef SIM3_USE_ROBUST_EDGES
	printf("reprojection error avg: %g px, RMSE: %g px (robustified, " PRIsize " observations, %.2f weight)\n",
		f_AVG_reprojection_err_robust, f_RMSE_reprojection_err_robust, edge_pool.n_Size(), f_robust_weight_sum);
#endif // SIM3_USE_ROBUST_EDGES
	printf("reprojection error avg: %g px, RMSE: %g px (visible points only, " PRIsize " observations)\n",
		f_AVG_reprojection_err_visible, f_RMSE_reprojection_err_visible, n_visible_obs_num);
	printf("there were " PRIsize " points fully visible, " PRIsize
		" points partially visible and " PRIsize " points fully invisible\n",
		n_fully_visible, n_partly_invisible, n_fully_invisible);
	// print also the visibility info

	//m_p_optimizer->r_Solver().Dump_ReprojectionErrors(); // debug

	m_p_optimizer->r_Solver().Dump();
	double f_error = m_p_optimizer->r_Solver().f_Chi_Squared_Error_Denorm();
	printf("denormalized chi2 error: %.2f\n", f_error);
	
	//m_p_optimizer->r_Solver().Dump_ReprojectionErrors(); // debug

	if(!m_p_optimizer->r_Solver().r_MarginalCovariance().r_SparseMatrix().b_Empty())
		m_p_optimizer->r_Solver().r_MarginalCovariance().Dump_Diagonal();

	{
		Eigen::VectorXd diag;
		m_p_optimizer->r_Solver().r_Lambda().Get_Diagonal(diag, true);
		double f_max_diag = diag.array().abs().maxCoeff();
		double f_min_diag = diag.array().abs().minCoeff();
		printf("approximate condition number of the information matrix is %.10g (min diag %g, max diag %g)\n",
			f_max_diag / f_min_diag, f_min_diag, f_max_diag);
	}
}

bool SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Dump_State(const char *p_s_filename) const
{
	return m_p_optimizer->r_System().Dump(p_s_filename);
}

/**
 *	@brief saves SE(3) state file from the optimizer system
 */
bool SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Dump_Marginals(const char *p_s_filename) const
{
	if(!m_p_optimizer->r_Solver().r_MarginalCovariance().r_SparseMatrix().b_Empty())
		return m_p_optimizer->r_Solver().r_MarginalCovariance().Dump_Diagonal(p_s_filename);
	return true;
}

bool SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Dump_State_SE3(const char *p_s_filename) const
{
	const CBAOptimizerCore_Sim3_gXYZ_gXYZ::CSystemType::_TyVertexMultiPool &vertex_pool = m_p_optimizer->r_System().r_Vertex_Pool();

	FILE *p_fw;
	if(!(p_fw = fopen(p_s_filename, "w")))
		return false;

	CSE3StatePrimitiveCallback<CBAOptimizerCore_Sim3_gXYZ_gXYZ::_TyLandmark> state_writer(p_fw, *this);
	vertex_pool.For_Each(state_writer);

	if(ferror(p_fw)) {
		fclose(p_fw);
		return false;
	}

	fclose(p_fw);
	return true;
}

bool SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Dump_Poses_SE3(const char *p_s_filename) const
{
	const CBAOptimizerCore_Sim3_gXYZ_gXYZ::CSystemType::_TyVertexMultiPool &vertex_pool = m_p_optimizer->r_System().r_Vertex_Pool();

	FILE *p_fw;
	if(!(p_fw = fopen(p_s_filename, "w")))
		return false;

	CSE3StatePrimitiveCallback<CBAOptimizerCore_Sim3_gXYZ_gXYZ::_TyLandmark> state_writer(p_fw, *this, true);
	vertex_pool.For_Each(state_writer);

	if(ferror(p_fw)) {
		fclose(p_fw);
		return false;
	}

	fclose(p_fw);
	return true;
}

/**
 *	@brief saves a graph file from the optimizer system
 */
bool SlamPP_Optimizer_Sim3_gXYZ_gXYZ::Dump_Graph_SE3(const char *p_s_filename) const
{
	const CBAOptimizerCore_Sim3_gXYZ_gXYZ::CSystemType::_TyVertexMultiPool &vertex_pool = m_p_optimizer->r_System().r_Vertex_Pool();
	const CBAOptimizerCore_Sim3_gXYZ_gXYZ::CSystemType::_TyEdgeMultiPool &edge_pool = m_p_optimizer->r_System().r_Edge_Pool();

	FILE *p_fw;
	if(!(p_fw = fopen(p_s_filename, "w")))
		return false;

	CGraphPrimitiveCallback<CBAOptimizerCore_Sim3_gXYZ_gXYZ::_TyLandmark,CBAOptimizerCore_Sim3_gXYZ_gXYZ::_TyObservation> graph_writer(p_fw, *this);
	vertex_pool.For_Each(graph_writer);
	edge_pool.For_Each(graph_writer);

	if(ferror(p_fw)) {
		fclose(p_fw);
		return false;
	}

	fclose(p_fw);
	return true;
}

/*
 *								=== ~SlamPP_Optimizer_Sim3_gXYZ_gXYZ ===
 */
}
}
