#pragma once

#include <openMVG/vsslam/optimization/slampp/optimizers/SlamPP_Optimizer_Sim3_gXYZ_gXYZ.hpp>

namespace openMVG {
namespace vsslam {

class CCalcRepErr { // todo - maybe change this so that it can collect all the stats in a single For_Each() call
protected:
  double &m_r_f_dest;
  double &m_r_f_robust_weight;
  bool &m_r_b_visible;

public:
  CCalcRepErr(double &r_f_dest, double &r_f_robust_weight, bool &r_b_visible)
    :m_r_f_dest(r_f_dest), m_r_f_robust_weight(r_f_robust_weight), m_r_b_visible(r_b_visible)
  {}

  template <class _TyEdge>
  void operator ()(const _TyEdge &r_t_edge)
  {
    m_r_f_dest = r_t_edge.f_Reprojection_Error(m_r_b_visible);
#ifdef SIM3_USE_ROBUST_EDGES
    Eigen::Vector1d r; r(0) = m_r_f_dest; // can't initialize directly, the constructor is ambiguous
    m_r_f_robust_weight = r_t_edge.f_RobustWeight(r);
#else // SIM3_USE_ROBUST_EDGES
    m_r_f_robust_weight = 1;
#endif // SIM3_USE_ROBUST_EDGES
    // this function is not exported via the facade, needs to be called like this
  }

  void operator ()(const CEdgePoseCamSim3 &UNUSED(r_t_edge)) // ignore pose-pose edges for calculation of reprojection error
  {
    m_r_f_dest = 0;
    m_r_f_robust_weight = 0; // !!
  }
};

}
}
