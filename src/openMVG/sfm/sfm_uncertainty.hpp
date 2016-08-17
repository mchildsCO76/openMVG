// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_UNCERTAINTY_HPP
#define OPENMVG_SFM_UNCERTAINTY_HPP

#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/numeric/numeric.h"
#include "ceres/types.h"
#include "ceres/cost_function.h"

#include "openMVG/types.hpp"
#include <cereal/cereal.hpp> // Serialization

namespace openMVG {
namespace sfm {

struct SfM_Data;

struct ObservationUncertainty
{
  ObservationUncertainty():id_observation(UndefinedIndexT), covariance(Eigen::Matrix<double, 2, 2>::Zero()) { }
  ObservationUncertainty(IndexT idObservation): id_observation(idObservation), covariance(Eigen::Matrix<double, 2, 2>::Zero()) {}

  IndexT id_observation;
  Eigen::Matrix<double, 2, 2> covariance;

  // Serialization
  template <class Archive>
  void save( Archive & ar) const
  {
    ar(cereal::make_nvp("id_observation", id_observation ));

    const std::vector<std::vector<double>> mat =
    {
      { covariance( 0, 0 ), covariance( 0, 1 ),
       covariance( 1, 0 ), covariance( 1, 1 )}
    };
    ar( cereal::make_nvp( "covariance", mat ) );
  }

  // Serialization
  template <class Archive>
  void load( Archive & ar)
  {
    ar(cereal::make_nvp("id_observation", id_observation ));

    std::vector<std::vector<double>> mat( 2, std::vector<double>( 2 ) );
    ar( cereal::make_nvp( "covariance", mat ) );
    // copy back to the covariance
      covariance.row( 0 ) = Eigen::Map<const Eigen::Matrix<double, 2, 1>>( &( mat[0][0] ) );
      covariance.row( 1 ) = Eigen::Map<const Eigen::Matrix<double, 2, 1>>( &( mat[0][3] ) );
  }
};

struct PoseUncertainty
{
  PoseUncertainty():id_pose(UndefinedIndexT),covariance(Eigen::Matrix<double, 6, 6>::Zero()) { }
  PoseUncertainty(IndexT idPose): id_pose(idPose), covariance(Eigen::Matrix<double, 6, 6>::Zero()) {}

  IndexT id_pose;
  Eigen::Matrix<double, 6, 6> covariance;

  // Serialization
  template <class Archive>
  void save( Archive & ar) const
  {
    ar(cereal::make_nvp("id_pose", id_pose ));

    const std::vector<std::vector<double>> mat =
    {
      { covariance( 0, 0 ), covariance( 0, 1 ), covariance( 0, 2 ), covariance( 0, 3 ), covariance( 0, 4 ), covariance( 0, 5 ),
       covariance( 1, 0 ), covariance( 1, 1 ), covariance( 1, 2 ), covariance( 1, 3 ), covariance( 1, 4 ), covariance( 1, 5 ) ,
       covariance( 2, 0 ), covariance( 2, 1 ), covariance( 2, 2 ), covariance( 2, 3 ), covariance( 2, 4 ), covariance( 2, 5 ) ,
       covariance( 3, 0 ), covariance( 3, 1 ), covariance( 3, 2 ), covariance( 3, 3 ), covariance( 3, 4 ), covariance( 3, 5 ) ,
       covariance( 4, 0 ), covariance( 4, 1 ), covariance( 4, 2 ), covariance( 4, 3 ), covariance( 4, 4 ), covariance( 4, 5 ) ,
       covariance( 5, 0 ), covariance( 5, 1 ), covariance( 5, 2 ), covariance( 5, 3 ), covariance( 5, 4 ), covariance( 5, 5 ) }
    };
    ar( cereal::make_nvp( "covariance", mat ) );
  }

  // Serialization
  template <class Archive>
  void load( Archive & ar)
  {
    ar(cereal::make_nvp("id_pose", id_pose ));

    std::vector<std::vector<double>> mat( 6, std::vector<double>( 6 ) );
    ar( cereal::make_nvp( "covariance", mat ) );
    // copy back to the covariance
    covariance.row( 0 ) = Eigen::Map<const Eigen::Matrix<double, 6, 1> >( &( mat[0][0] ) );
    covariance.row( 1 ) = Eigen::Map<const Eigen::Matrix<double, 6, 1> >( &( mat[0][6] ) );
    covariance.row( 2 ) = Eigen::Map<const Eigen::Matrix<double, 6, 1> >( &( mat[0][12] ) );
    covariance.row( 3 ) = Eigen::Map<const Eigen::Matrix<double, 6, 1> >( &( mat[0][18] ) );
    covariance.row( 4 ) = Eigen::Map<const Eigen::Matrix<double, 6, 1> >( &( mat[0][24] ) );
    covariance.row( 5 ) = Eigen::Map<const Eigen::Matrix<double, 6, 1> >( &( mat[0][30] ) );
  }
};

struct LandmarkUncertainty
{
  LandmarkUncertainty():id_landmark(UndefinedIndexT),covariance(Eigen::Matrix<double, 3, 3>::Zero()) { }
  LandmarkUncertainty(IndexT idLandmark): id_landmark(idLandmark), covariance(Eigen::Matrix<double, 3, 3>::Zero()) {}

  IndexT id_landmark;
  Eigen::Matrix<double, 3, 3> covariance;

  // Serialization
  template <class Archive>
  void save( Archive & ar) const
  {
    ar(cereal::make_nvp("id_landmark", id_landmark ));

    const std::vector<std::vector<double>> mat =
    {
      { covariance( 0, 0 ), covariance( 0, 1 ), covariance( 0, 2 ) ,
       covariance( 1, 0 ), covariance( 1, 1 ), covariance( 1, 2 ) ,
       covariance( 2, 0 ), covariance( 2, 1 ), covariance( 2, 2 ) }
    };
    ar( cereal::make_nvp( "covariance", mat ) );
  }

  // Serialization
  template <class Archive>
  void load( Archive & ar)
  {
    ar(cereal::make_nvp("id_landmark", id_landmark ));

    std::vector<std::vector<double>> mat( 3, std::vector<double>( 3 ) );
    ar( cereal::make_nvp( "covariance", mat ) );
    // copy back to the covariance
      covariance.row( 0 ) = Eigen::Map<const Vec3>( &( mat[0][0] ) );
      covariance.row( 1 ) = Eigen::Map<const Vec3>( &( mat[0][3] ) );
      covariance.row( 2 ) = Eigen::Map<const Vec3>( &( mat[0][6] ) );
  }
};

bool EstimateUncertainty
(
  // the SfM scene to refine
  SfM_Data & sfm_data,
  // Bundle adjustment object
  Bundle_Adjustment_Ceres &bundle_adjustment_obj,
  // tell which parameter needs to be adjusted
  const Optimize_Options &options,
  const bool estimateLandmarks = false
);

void EstimateQualityOfStructure
(
  const SfM_Data & sfm_data,
  std::vector<double> & vec_structureUncertainty
);

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_UNCERTAINTY_HPP
