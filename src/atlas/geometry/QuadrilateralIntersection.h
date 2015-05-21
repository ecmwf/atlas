/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_QuadIntersection_h
#define atlas_QuadIntersection_h

#include <limits>

#include "eckit/eckit_config.h"

#ifdef HAVE_EIGEN

#include "eckit/maths/Eigen.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/geometry/Ray.h"

namespace atlas {
namespace geometry {

//----------------------------------------------------------------------------------------------------------------------

class QuadrilateralIntersection {

public:

  QuadrilateralIntersection(const double* x0, const double* x1, const double* x2, const double* x3) {

    v00 = Eigen::Vector3d::Map(x0);
    v10 = Eigen::Vector3d::Map(x1);
    v11 = Eigen::Vector3d::Map(x2);
    v01 = Eigen::Vector3d::Map(x3);

  }

  Intersect intersects(const Ray& r, double epsilon = 5 * std::numeric_limits<double>::epsilon()) const;

private: // members

  Eigen::Vector3d v00;
  Eigen::Vector3d v10;
  Eigen::Vector3d v11;
  Eigen::Vector3d v01;

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif  // HAVE_EIGEN

#endif
