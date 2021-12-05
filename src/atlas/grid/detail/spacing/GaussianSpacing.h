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

#include <string>

#include "atlas/grid/detail/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace spacing {

// clang-format off

/// @brief Gaussian spacing in interval
///
/// There are N Gaussian spaced points in the open interval (90, -90)
///
/// Using the constructor GaussianSpacing( N ) we can create
///
///    GaussianSpacing( 4 )                -->  { 59.44... , 19.87... , -19.87... , -59.44... }
///
/// Configuration parameters can be passed as well with following keys:
///
///    {"N":4 }                            -->  { 59.44... , 19.87... , -19.87... , -59.44... }
///
/// To reverse the orientation of points to go from negative to positive
/// instead, pass also the start and end keys:
///
///    {"N":4, "start":-90, "end":90 }    -->  { -59.44... , -19.87... ,  19.87... , 59.44... }

// clang-format on

class GaussianSpacing : public Spacing {
public:
    // constructor
    GaussianSpacing(const eckit::Parametrisation& p);

    GaussianSpacing(long N);

    // class name
    static std::string static_type() { return "gaussian"; }
    virtual std::string type() const { return static_type(); }

    virtual Spec spec() const;
};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
