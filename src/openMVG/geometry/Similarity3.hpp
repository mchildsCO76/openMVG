// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GEOMETRY_SIMILARITY3_HPP
#define OPENMVG_GEOMETRY_SIMILARITY3_HPP

#include "openMVG/geometry/pose3.hpp"

namespace openMVG
{
namespace geometry
{

/**
* @brief Define a 3D Similarity transform encoded as a 3D pose plus a scale
*/
struct Similarity3
{
  /// Pose
  Pose3 pose_;

  /// Scale
  double scale_;

  /**
  * @brief Default constructor
  * @note This define identity transformation centered at origin of a cartesian frame
  */
  Similarity3()
    : pose_( Pose3() ),
      scale_( 1.0 )
  {

  }

  /**
  * @brief Constructor
  * @param pose a 3d pose
  * @param scale a scale factor
  */
  Similarity3( const Pose3 & pose, const double scale )
    : pose_( pose ),
      scale_( scale )
  {

  }

  /**
  * @brief Get Rotation matrix
  * @return Rotation matrix
  */
  const Mat3& rotation() const
  {
    return pose_.rotation();
  }

  /**
  * @brief Get Rotation matrix
  * @return Rotation matrix
  */
  Mat3& rotation()
  {
    return pose_.rotation();
  }

  /**
  * @brief Get center of rotation
  * @return center of rotation
  */
  const Vec3& center() const
  {
    return pose_.center();
  }

  /**
  * @brief Get center of rotation
  * @return Center of rotation
  */
  Vec3& center()
  {
    return pose_.center();
  }

  /**
  * @brief Get translation vector
  * @return translation vector
  * @note t = -RC
  */
  inline Vec3 translation() const
  {
    return pose_.translation();
  }

  double& scale()
  {
    return scale_;
  }

  const double& scale() const
  {
    return scale_;
  }

  /**
  * @brief Apply transformation to a point
  * @param point Input point
  * @return transformed point
  */
  Mat3X operator () ( const Mat3X & point ) const
  {
    return scale_ * pose_( point );
  }

  /**
  * @brief Concatenation of pose
  * @param pose Pose to be concatenated with the current one
  * @return Concatenation of poses
  */
  Pose3 operator () ( const Pose3 & pose ) const
  {
    return Pose3( pose.rotation() * pose_.rotation().transpose(), this->operator()( pose.center() ) );
  }

  /**
  * @brief Get inverse of the similarity
  * @return Inverse of the similarity
  */
  Similarity3 inverse() const
  {
    return Similarity3(pose_.inverse(), 1.0 / scale_);
  }



  Mat4 transformation() const
  {
    Mat4 T;
    T.block(0,0,3,3) << scale_*pose_.rotation();
    T.block(0,3,3,1) << pose_.translation();
    T(3,3) = 1;
    return T;
  }



};

} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_SIMILARITY3_HPP
