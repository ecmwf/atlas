/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/eckit_config.h"
#include "eckit/exception/Exceptions.h"

#ifdef HAVE_EIGEN

#include "eckit/maths/Eigen.h"

#include "atlas/QuadIntersection.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

bool Quad::intersects(const Ray& r, Isect& isect, const double slack) {

  NOTIMP;

  return true;
}

//------------------------------------------------------------------------------------------------------

}  // namespace atlas

#endif  // HAVE_EIGEN
