#pragma once

#define SIM3_USE_ROBUST_EDGES
#define USE_EXPLICIT_HARD_UF // use explicit UF via a const vertex?
#include <string.h>
#include <stdio.h>


#include "slam/LinearSolver_UberBlock.h"
//#include "slam/LinearSolver_CholMod.h" // not used
#include "slam/ConfigSolvers.h" // nonlinear graph solvers
//#include "slam/Timer.h" // not used
#include "slam/3DSolverBase.h" // want C3DJacobians::Quat_to_AxisAngle() and C3DJacobians::AxisAngle_to_Quat()
#include "slam/Sim3SolverBase.h" // C3DJacobians, CBAJacobians, generally useful functions for BA and SE(3), does not need to be included

#include "slam/Sim3_Types.h"
//#include "slam/Marginals.h" // included from nonlinear solver, if supported
#include "slam/NonlinearSolver_Lambda_DL.h"
//#include "slam/BASolverBase.h" // included from BA_Types.h
#include "slam/Eigenvalues.h"


namespace openMVG {
namespace VSSLAM {


class CBAOptimizerCore_Sim3_gXYZ_gXYZ {
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

  #ifdef USE_EXPLICIT_HARD_UF
    using TConstVertexTypelist = MakeTypelist_Safe((CVertexCamSim3));
    using TEdgeTypelist = MakeTypelist_Safe((_TyObservation, CEdgePoseCamSim3));
    using CSystemType = CFlatSystem<CBaseVertex, TVertexTypelist, CBaseEdge, TEdgeTypelist, CNullUnaryFactorFactory,
        CVertexCamSim3, TConstVertexTypelist>;
  #else // USE_EXPLICIT_HARD_UF
    using TEdgeTypelist = MakeTypelist_Safe((_TyObservation));
    using CSystemType = CFlatSystem<CBaseVertex, TVertexTypelist,TEdgeTypelist>;
  #endif // USE_EXPLICIT_HARD_UF


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


}
}
