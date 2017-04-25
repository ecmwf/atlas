/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_interpolation_element_Quad2D_h
#define atlas_interpolation_element_Quad2D_h

#include <limits>

#include "atlas/interpolation/Vector2D.h"
#include "atlas/interpolation/method/Intersect.h"

namespace atlas {
namespace interpolation {
namespace element {

//----------------------------------------------------------------------------------------------------------------------

class Quad2D {
public:

    Quad2D(const double* x0, const double* x1, const double* x2, const double* x3) {
        v00 = Vector2D::Map(x0);
        v10 = Vector2D::Map(x1);
        v11 = Vector2D::Map(x2);
        v01 = Vector2D::Map(x3);
    }

    method::Intersect intersects(
            const Vector2D& p,
            double edgeEpsilon = 5 * std::numeric_limits<double>::epsilon(),
            double epsilon = 5 * std::numeric_limits<double>::epsilon() ) const;

    bool validate() const;

    double area() const;

    void print(std::ostream& s) const {
        s << "Quad2D[v00=" << v00 << ",v10=" << v10 << ",v11=" << v11 << ",v01=" << v01 << "]";
    }

    friend std::ostream& operator<<(std::ostream& s, const Quad2D& p) {
        p.print(s);
        return s;
    }

private:  // members
    Vector2D v00; // aka v0
    Vector2D v10; // aka v1
    Vector2D v11; // aka v2
    Vector2D v01; // aka v3

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas


#endif
