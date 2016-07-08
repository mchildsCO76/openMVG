// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_UNCERTAINTY_HPP
#define OPENMVG_SFM_UNCERTAINTY_HPP

#include "openMVG/types.hpp"
#include "openMVG/numeric/numeric.h"

#include <cereal/cereal.hpp> // Serialization

namespace openMVG {
namespace sfm {


/// Define a landmark (a 3D point, with it's 2d observations)
struct UncertaintyLandmark
{
  Eigen::Matrix3d covariance;
  double meanReprojError;

  // Serialization
  template <class Archive>
  void save( Archive & ar) const
  {
    const std::vector<double> cov = { covariance(0,0), covariance(0,1), covariance(0,2), covariance(1,0), covariance(1,1), covariance(1,2), covariance(2,0), covariance(2,1), covariance(2,2) };
    ar(cereal::make_nvp("covariance", cov ));
    ar(cereal::make_nvp("meanReprojError", meanReprojError ));
  }

  template <class Archive>
  void load( Archive & ar)
  {
    std::vector<double> cov(9);
    ar(cereal::make_nvp("covariance", cov ));
    covariance = Eigen::Map<const Eigen::Matrix3d>(&cov[0]);
    ar(cereal::make_nvp("meanReprojError", meanReprojError ));
  }
};

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> EigenSparseMatrix;

struct UncertaintyCams
{
  EigenSparseMatrix covariance;

  // Serialization
  template <class Archive>
  void save( Archive & ar) const
  {
	const int n_rows = covariance.rows();
	const int n_cols = covariance.cols();
	ar(cereal::make_nvp("rows", n_rows ));
	ar(cereal::make_nvp("cols", n_cols ));
	
	// CRS Matrix
	std::vector<int> m_rows;
	std::vector<int> m_cols;
	std::vector<double>m_values;
	
	int n_element = 0;
	for(int i=0; i < n_rows; ++i) {
	  m_rows.push_back(n_element);
	  for(typename EigenSparseMatrix::InnerIterator it(covariance,i); it; ++it) {
	    m_cols.push_back(it.col());
	    m_values.push_back(it.value());
	    n_element++;
	  }
	}
	
	const std::vector<int> r_const(m_rows.begin(),m_rows.end());
  const std::vector<int> c_const(m_cols.begin(),m_cols.end());
  const std::vector<double> v_const(m_values.begin(),m_values.end());

	ar(cereal::make_nvp("n_elements", n_element ));
	ar(cereal::make_nvp("rows_data", r_const ));
	ar(cereal::make_nvp("cols_data", c_const ));
	ar(cereal::make_nvp("values_data", v_const ));
	
  }

  template <class Archive>
  void load( Archive & ar)
  {
	int n_rows;
	int n_cols;
	int n_element;
	ar(cereal::make_nvp("rows", n_rows ));
	ar(cereal::make_nvp("cols", n_cols ));
	ar(cereal::make_nvp("n_elements", n_element ));
	
	std::vector<int> m_rows(n_rows+1);
	std::vector<int> m_cols(n_element);
	std::vector<double>m_values(n_element);
  ar(cereal::make_nvp("rows_data", m_rows ));
  ar(cereal::make_nvp("cols_data", m_cols ));
  ar(cereal::make_nvp("values_data", m_values ));

  covariance = Eigen::MappedSparseMatrix<double, Eigen::RowMajor>(
  	n_rows, n_cols,
  	static_cast<int>(n_element),
  	m_rows.data(), m_cols.data(), m_values.data());

  }
};



} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_UNCERTAINTY_HPP
