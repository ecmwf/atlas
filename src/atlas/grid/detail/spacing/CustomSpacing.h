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

#include <array>
#include <string>

#include "atlas/grid/detail/spacing/Spacing.h"

namespace atlas {
namespace grid {
namespace spacing {


// clang-format off

/// @brief Custom spacing in interval
///
/// There are N points, by default in the open interval (90, -90). The interval is used to
/// determine if a grid's domain contains poles.
///
/// Using the constructor CustomSpacing( N, x[], {min,max} ) we can create
///
///    CustomSpacing( 4, {75,25,-25,-75} )                  -->  { 75 , 25 , -25, -75 }
///    CustomSpacing( 4, {75,25,-25,-75}, {-90,90} )        -->  { 75 , 25 , -25, -75 }
///
/// The optional argument {min,max} serves as purpose to indicate that the points
/// lie in the open interval (min,max). If not specified, the default values are taken
/// to be the North and South pole's latitudes.
///
/// Configuration parameters can be passed as well with following keys:
///
///    {"N":4, "values":[75,25,-25,75] }                        -->  { 75 , 25 , -25 , -75 }
///    {"N":4, "values":[75,25,-25,75], "interval":[-90,90] }   -->  { 75 , 25 , -25 , -75 }

// clang-format on

class CustomSpacing : public Spacing {
private:
    using Interval = std::array<double, 2>;
    static constexpr double North() { return 90.; }
    static constexpr double South() { return -90.; }

public:
    // constructor
    CustomSpacing(const eckit::Parametrisation& p);

    CustomSpacing(long N, const double x[], const Interval& = {North(), South()});

    // class name
    static std::string static_type() { return "custom"; }
    virtual std::string type() const { return static_type(); }

    virtual Spec spec() const;
};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
