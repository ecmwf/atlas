/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/geometry/QuadrilateralIntersection.h"

#include <cmath>

#include "eckit/eckit_config.h"
#include "eckit/exception/Exceptions.h"

#ifdef HAVE_EIGEN

#include "eckit/maths/Eigen.h"
#include "eckit/log/Log.h"

#include "atlas/geometry/TriangleIntersection.h"

using eckit::Log;
using Eigen::Vector3d;

namespace atlas {
namespace geometry {

//----------------------------------------------------------------------------------------------------------------------


Intersect QuadrilateralIntersection::intersects(const Ray &r, double epsilon) const
{
    Intersect isect; // intersection is false

    TriangleIntersection T013(v00.data(), v10.data(), v01.data());
    isect = T013.intersects(r,epsilon);
    if(isect)
        return isect;

    TriangleIntersection T231(v11.data(), v01.data(), v10.data());
    isect = T231.intersects(r,epsilon);
    if(isect)
    {
        isect.u = 1 - isect.u;
        isect.v = 1 - isect.v;
        return isect;
    }

    return isect.success(false);
}

bool QuadrilateralIntersection::validate() const {

    // normal for sub-triangle T231

    Vector3d E23 = v01 - v11;
    Vector3d E21 = v10 - v11;

    Vector3d N231 = E23.cross(E21);

    // normal for sub-triangle T013

    Vector3d E01 = v10 - v00;
    Vector3d E03 = v01 - v00;

    Vector3d N013 = E01.cross(E03);

    // normal for sub-triangle T120

    Vector3d E12 = - E21;
    Vector3d E10 = - E01;

    Vector3d N120 = E12.cross(E10);

    // normal for sub-triangle T302

    Vector3d E30 = - E03;
    Vector3d E32 = - E23;

    Vector3d N302 = E30.cross(E32);

    // all normals must point same way

    double dot02 = N231.dot(N013);
    double dot23 = N013.dot(N120);
    double dot31 = N120.dot(N302);
    double dot10 = N302.dot(N231);

    // all normals must point same way

    if (( dot02 >= 0. &&  dot23 >= 0. &&  dot31 >= 0. &&  dot10 >= 0. ) ||
        ( dot02 <= 0. &&  dot23 <= 0. &&  dot31 <= 0. &&  dot10 <= 0. ))
        return true;
    else
        return false;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace geometry
}  // namespace atlas

#endif  // HAVE_EIGEN
