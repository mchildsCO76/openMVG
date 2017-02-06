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
  const Vec3 center() const
  {
    return pose_.center();
  }

  /**
  * @brief Get center of rotation
  * @return Center of rotation
  */
  Vec3 center()
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
  * @brief Concatenation of pose
  * @param pose Pose to be concatenated with the current one
  * @return Concatenation of poses
  */
  Similarity3 operator () ( const Similarity3 & pose ) const
  {
    return Similarity3( Pose3(pose_.rotation()* pose.pose_.rotation(), pose.pose_.center() + (1/pose.scale_) * pose.pose_.rotation().transpose() * pose_.center()), scale_*pose.scale_ );
  }

  /**
  * @brief Get inverse of the similarity
  * @return Inverse of the similarity
  */
  Similarity3 inverse() const
  {
    return Similarity3(pose_.inverse(), 1.0 / scale_);
  }


  /**
  * @brief Return the depth (distance) of a point respect to the camera center
  * @param X Input point
  * @return Distance to center
  */
  double depth( const Vec3 &X ) const
  {
    return ( scale_ * pose_.rotation() * ( X - pose_.center() ) )[2];
  }

};

} // namespace geometry
} // namespace openMVG

#endif  // OPENMVG_GEOMETRY_SIMILARITY3_HPP
