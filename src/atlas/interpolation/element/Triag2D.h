/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_interpolation_element_Triag2D_h
#define atlas_interpolation_element_Triag2D_h

#include <limits>

#include "atlas/interpolation/Vector2D.h"
#include "atlas/interpolation/method/Intersect.h"

namespace atlas {
namespace interpolation {
namespace method {
struct Ray;
}
namespace element {

//----------------------------------------------------------------------------------------------------------------------

class Triag2D {

public: // types

    Triag2D(const Vector2D& x0, const Vector2D& x1, const Vector2D& x2):
        v0(x0),
        v1(x1),
        v2(x2) {
    }

    Triag2D(const double* x0, const double* x1, const double* x2) {
        v0 = Vector2D::Map(x0);
        v1 = Vector2D::Map(x1);
        v2 = Vector2D::Map(x2);
    }

    method::Intersect intersects(
            const method::Ray& r,
            double edgeEpsilon = 5 * std::numeric_limits<double>::epsilon(),
            double epsilon = 5 * std::numeric_limits<double>::epsilon() ) const;

    double area() const;

    void print(std::ostream& s) const {
        s << "Triag2D["
          <<  "v0=("   << v0[0] << ", " << v0[1]
          << "), v1=(" << v1[0] << ", " << v1[1]
          << "), v2=(" << v2[0] << ", " << v2[1] << ")]";
    }

    friend std::ostream& operator<<(std::ostream& s, const Triag2D& p) {
        p.print(s);
        return s;
    }

private: // members

    Vector2D v0;
    Vector2D v1;
    Vector2D v2;

};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace element
}  // namespace interpolation
}  // namespace atlas


#endif
