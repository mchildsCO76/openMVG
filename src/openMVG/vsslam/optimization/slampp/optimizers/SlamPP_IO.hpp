#pragma once

#include <openMVG/vsslam/optimization/slampp/optimizers/SlamPP_Optimizer_Sim3_gXYZ_gXYZ.hpp>

namespace openMVG {
namespace vsslam {

template <typename _TyLandmark, typename _TyObservation>
class CGraphPrimitiveCallback {
protected:
  size_t m_n_id; /**< @brief zero-based id of the last written vertex @note This is needed because the vertices do not store their id, but they will be enumerated in the order they go. */
  FILE *m_p_fw; /**< @brief output file */
  const SlamPP_Optimizer &m_r_optimizer; /**< @brief reference to the optimizer */

public:
  inline CGraphPrimitiveCallback(FILE *p_fw, const SlamPP_Optimizer &r_optimizer)
    :m_p_fw(p_fw), m_n_id(size_t(-1)), m_r_optimizer(r_optimizer)
  {}

  void operator ()(const _TyLandmark &r_vertex_invd)
  {
    Eigen::Vector3d v_state = m_r_optimizer.v_XYZVertex(++ m_n_id);
    fprintf(m_p_fw, "VERTEX_XYZ " PRIsize " %.15f %.15f %.15f\n", m_n_id, v_state(0), v_state(1), v_state(2));
  }

  void operator ()(const CVertexCamSim3 &r_vertex_cam)
  {
    const Eigen::Matrix<double, 7, 1> &r_v_state = r_vertex_cam.r_v_State();
    const Eigen::Matrix<double, 5, 1> &r_v_intrinsics = r_vertex_cam.v_Intrinsics();
    // get state and intrinsics

    CSim3Jacobians::TSim3 t_pose(r_v_state, CSim3Jacobians::TSim3::from_sim3_vector);
    Eigen::Quaterniond t_rotation = t_pose.t_rotation;
    Eigen::Vector3d v_translation = t_pose.v_translation / t_pose.f_scale; // fixme - is this correct?
    // convert to SE(3), mind the scale this time

    fprintf(m_p_fw, "VERTEX_CAM " PRIsize " %.15f %.15f %.15f %.15g %.15g %.15g %.15g %.15f %.15f %.15f %.15f %.15f\n", // use %g for the quat, as it will have less than zero numbers
      ++ m_n_id, v_translation(0), v_translation(1), v_translation(2), t_rotation.x(), t_rotation.y(),
      t_rotation.z(), t_rotation.w(), r_v_intrinsics(0), r_v_intrinsics(1), r_v_intrinsics(2),
      r_v_intrinsics(3), r_v_intrinsics(4));
  }

  void operator ()(const _TyObservation &r_t_edge)
  {
    size_t n_point_id = r_t_edge.n_Vertex_Id(0);
    size_t n_camera_id = r_t_edge.n_Vertex_Id(1); // not swapped

    const CEdgeP2CSim3G::_TyVectorAlign &r_v_measurement = r_t_edge.v_Measurement();
    const CEdgeP2CSim3G::_TyMatrixAlign &r_t_information = r_t_edge.t_Sigma_Inv();

    fprintf(m_p_fw, "EDGE_PROJECT_P2MC " PRIsize " " PRIsize " %f %f %g %g %g\n", // these do not need to be very precise (x, y detected by SIFT with barely any fractional precision, covariance typically large integer numbers)
      n_point_id, n_camera_id, r_v_measurement(0), r_v_measurement(1),
      r_t_information(0, 0), r_t_information(0, 1), r_t_information(1, 1)); // only the upper triangular part is stored
  }

  void operator ()(const CEdgePoseCamSim3 &r_t_edge)
  {
    size_t n_point0_id = r_t_edge.n_Vertex_Id(0);
    size_t n_point1_id = r_t_edge.n_Vertex_Id(1);

    fprintf(m_p_fw, "# hard UF edge " PRIsize " " PRIsize "\n", // these do not need to be very precise (x, y detected by SIFT with barely any fractional precision, covariance typically large integer numbers)
      n_point0_id, n_point1_id); // only the upper triangular part is stored
  }
};

template <typename _TyLandmark>
class CSE3StatePrimitiveCallback {
protected:
  size_t m_n_id; /**< @brief zero-based id of the last written vertex @note This is needed because the vertices do not store their id, but they will be enumerated in the order they go. */
  FILE *m_p_fw; /**< @brief output file */
  const SlamPP_Optimizer &m_r_optimizer; /**< @brief reference to the optimizer */
  bool m_b_ignore_landmarks; /**< @brief landmarks ignore flag */

public:
  inline CSE3StatePrimitiveCallback(FILE *p_fw, const SlamPP_Optimizer &r_optimizer, bool b_ignore_landmarks = false)
    :m_p_fw(p_fw), m_n_id(size_t(-1)), m_r_optimizer(r_optimizer), m_b_ignore_landmarks(b_ignore_landmarks)
  {}

  void operator ()(const _TyLandmark &r_vertex_invd)
  {
    if(m_b_ignore_landmarks)
      return;
    Eigen::Vector3d v_state = m_r_optimizer.v_XYZVertex(++ m_n_id);
    fprintf(m_p_fw, "%.15f %.15f %.15f\n", v_state(0), v_state(1), v_state(2));
  }

  void operator ()(const CVertexCamSim3 &r_vertex_cam)
  {
    ++ m_n_id; // !!

    const Eigen::Matrix<double, 7, 1> &r_v_state = r_vertex_cam.r_v_State();
    const Eigen::Matrix<double, 5, 1> &r_v_intrinsics = r_vertex_cam.v_Intrinsics();
    // get state and intrinsics

    CSim3Jacobians::TSim3 t_pose(r_v_state, CSim3Jacobians::TSim3::from_sim3_vector);
    t_pose.Invert(); // the SE(3) internal representation uses inverse camera poses. these do not get saved as inverse to the graph file, but they do get dumped like that in the internal format
    /*Eigen::Quaterniond t_rotation = t_pose.t_rotation;
    Eigen::Vector3d v_rotation;
    C3DJacobians::Quat_to_AxisAngle(t_rotation, v_rotation);
    Eigen::Vector3d v_translation = t_pose.v_translation / t_pose.f_scale;*/ // fixme - is this correct? // division seems to be correct
    Eigen::Vector7d v_trs = t_pose.v_tRs();
    v_trs.head<3>() /= v_trs(6); // fixme - is this correct? // division seems to be correct
    // convert to SE(3), mind the scale this time

    fprintf(m_p_fw, "%.15f %.15f %.15f %.15g %.15g %.15g\n", // use %g for the rotation, as it will have less than pi numbers
    //  v_translation(0), v_translation(1), v_translation(2), v_rotation(0), v_rotation(1), v_rotation(2));
      v_trs(0), v_trs(1), v_trs(2), v_trs(3), v_trs(4), v_trs(5));
  }
};

}
}
