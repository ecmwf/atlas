/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <iosfwd>

namespace atlas {
namespace interpolation {
namespace method {

//----------------------------------------------------------------------------------------------------------------------

/// Intersection data structure

struct Intersect {
    double u;
    double v;
    double t;

    Intersect();

    operator bool() const { return success_; }

    Intersect& success() {
        success_ = true;
        return *this;
    }
    Intersect& fail() {
        success_ = false;
        return *this;
    }

    void print(std::ostream& s) const;

    friend std::ostream& operator<<(std::ostream& s, const Intersect& p) {
        p.print(s);
        return s;
    }

private:
    bool success_;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
