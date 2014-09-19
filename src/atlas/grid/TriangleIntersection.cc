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

#ifndef HAVE_EIGEN
void atlas_grid_TriangeIntersection() {} // dummy symbol
#else

#include "eckit/maths/Eigen.h"

#include "atlas/grid/TriangleIntersection.h"

using namespace std;
using namespace atlas;

//------------------------------------------------------------------------------------------------------

namespace atlas {

bool triag_intersection( const Triag& tg, const Ray& r, Isect& isect, const double slack )
{
    Eigen::Vector3d edge1 = tg.v1 - tg.v0;
    Eigen::Vector3d edge2 = tg.v2 - tg.v0;
    Eigen::Vector3d pvec  = r.dir.cross(edge2);
    const double det = edge1.dot(pvec);
    if( std::abs(det) < slack ) return false;
    const double invDet = 1. / det;
    Eigen::Vector3d tvec = r.orig - tg.v0;
    isect.u = tvec.dot(pvec) * invDet;
    if( isect.u + slack < 0 || isect.u - slack > 1) return false;
    Eigen::Vector3d qvec = tvec.cross(edge1);
    isect.v = r.dir.dot(qvec) * invDet;
    if (isect.v + slack < 0 || isect.u + isect.v - slack > 1) return false;
    isect.t = edge2.dot(qvec) * invDet;
    return true;
}

//------------------------------------------------------------------------------------------------------

} // namespace eckit

#endif
