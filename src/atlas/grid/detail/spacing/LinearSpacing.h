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

/// @brief Linear spacing in interval
///
/// There are N equally spaced points in the closed interval [start, stop] or the
/// half-open interval [start, stop) (depending on whether endpoint is True or False)
///
/// Using the constructor LinearSpacing( start, end, N, endpoint ) we can create
///    LinearSpacing( 2, 3, 5, true  )                    -->  { 2.0 , 2.25 , 2.5 , 2.75 , 3.0 }
///    LinearSpacing( 2, 3, 5, false )                    -->  { 2.0 , 2.2  , 2.4 , 2.6  , 2.8 }
///
/// Configuration parameters can be passed as well with following keys:
///
///    {"start":2 , "end":3, "N":5, "endpoint":true }     --> { 2.0 , 2.25 , 2.5 , 2.75 , 3.0 }
///    {"start":2 , "end":3, "N":5, "endpoint":false}     --> { 2.0 , 2.2  , 2.4 , 2.6  , 2.8 }
///
/// Instead of the "end" key, you can provide the "length" key, to achieve the same results:
///
///    {"start":2 , "length":1, "N":5, "endpoint":true }  --> { 2.0 , 2.25 , 2.5 , 2.75 , 3.0 }
///    {"start":2 , "length":1, "N":5, "endpoint":false}  --> { 2.0 , 2.2  , 2.4 , 2.6  , 2.8 }

// clang-format on

class LinearSpacing : public Spacing {
public:
    using Interval = std::array<double, 2>;

public:
    /// constructor
    LinearSpacing(const eckit::Parametrisation& p);
    LinearSpacing(double start, double end, long N, bool endpoint = true);
    LinearSpacing(const Interval&, long N, bool endpoint = true);

    // LinearSpacing( double centre, double step, long N, bool endpoint=true );

    // class name
    static std::string static_type() { return "linear"; }
    virtual std::string type() const { return static_type(); }

    double step() const;

    bool endpoint() const;

    virtual Spec spec() const;

public:
    struct Params {
        double start;
        double end;
        long N;
        double length;
        bool endpoint;
        double step;
        Params() = default;
        Params(const eckit::Parametrisation& p);
        Params(double start, double end, long N, bool endpoint);
    };

protected:
    // points are equally spaced between xmin and xmax
    // Depending on value of endpoint, the spacing will be different
    void setup(double start, double end, long N, bool endpoint);

    double start_;
    double end_;
    long N_;
    bool endpoint_;
};

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
