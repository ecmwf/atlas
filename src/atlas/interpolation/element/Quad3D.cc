/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <iostream>

#include "atlas/interpolation/element/Quad3D.h"
#include "atlas/interpolation/element/Triag3D.h"
#include "atlas/interpolation/method/Ray.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace interpolation {
namespace element {

//----------------------------------------------------------------------------------------------------------------------

method::Intersect Quad3D::intersects( const method::Ray& r, double edgeEpsilon, double epsilon ) const {
    method::Intersect isect;  // intersection is false

    Triag3D T013( v00, v10, v01 );
    isect = T013.intersects( r, edgeEpsilon, epsilon );
    if ( isect ) {
        return isect;
    }

    Triag3D T231( v11, v01, v10 );
    isect = T231.intersects( r, edgeEpsilon, epsilon );
    if ( isect ) {
        isect.u = 1 - isect.u;
        isect.v = 1 - isect.v;
        return isect;
    }

    return isect.fail();
}

bool Quad3D::validate() const {
    // normal for sub-triangle T231

    Vector3D E23 = v01 - v11;
    Vector3D E21 = v10 - v11;

    Vector3D N231 = E23.cross( E21 );

    // normal for sub-triangle T013

    Vector3D E01 = v10 - v00;
    Vector3D E03 = v01 - v00;

    Vector3D N013 = E01.cross( E03 );

    // normal for sub-triangle T120

    Vector3D E12 = -E21;
    Vector3D E10 = -E01;

    Vector3D N120 = E12.cross( E10 );

    // normal for sub-triangle T302

    Vector3D E30 = -E03;
    Vector3D E32 = -E23;

    Vector3D N302 = E30.cross( E32 );

    // all normals must point same way

    double dot02 = N231.dot( N013 );
    double dot23 = N013.dot( N120 );
    double dot31 = N120.dot( N302 );
    double dot10 = N302.dot( N231 );

    // all normals must point same way
    bool is_inside = ( ( dot02 >= 0. && dot23 >= 0. && dot31 >= 0. && dot10 >= 0. ) ||
                       ( dot02 <= 0. && dot23 <= 0. && dot31 <= 0. && dot10 <= 0. ) );
    return is_inside;
}

double Quad3D::area() const {
    Triag3D T013( v00, v10, v01 );
    Triag3D T231( v11, v01, v10 );

    return T013.area() + T231.area();
}

void Quad3D::print( std::ostream& s ) const {
    auto printVector3D = [&s]( const Vector3D& v ) { s << "[" << v[0] << "," << v[1] << "," << v[2] << "]"; };
    s << "Quad3D[";
    printVector3D( v00 );
    s << ", ";
    printVector3D( v10 );
    s << ", ";
    printVector3D( v11 );
    s << ", ";
    printVector3D( v01 );
    s << "]";
}


//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas
