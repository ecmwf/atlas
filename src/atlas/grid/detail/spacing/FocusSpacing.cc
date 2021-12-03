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

#include "eckit/config/Parametrisation.h"

#include "atlas/grid/detail/spacing/FocusSpacing.h"
#include "atlas/grid/detail/spacing/SpacingFactory.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {
namespace spacing {

FocusSpacing::FocusSpacing(const eckit::Parametrisation& params) {
    double xmin;
    double xmax;
    long N;

    // retrieve xmin, xmax and N from params
    if (!params.get("start", xmin)) {
        throw_Exception("start missing in Params", Here());
    }
    if (!params.get("end", xmax)) {
        throw_Exception("end missing in Params", Here());
    }
    if (!params.get("N", N)) {
        throw_Exception("N missing in Params", Here());
    }

    // additional parameters for focus spacing
    if (!params.get("focus_factor", focus_factor_)) {
        throw_Exception("focus_factor missing in Params", Here());
    }

    x_.resize(N);
    if (N == 1) {
        x_[0] = 0.5 * (xmin + xmax);
    }
    else {
        const double midpoint = 0.5 * (xmin + xmax);
        const double d2       = 2. / double(N - 1);
        const double c1       = (xmax - xmin) * M_1_PI;
        const double c2       = 1. / focus_factor_;
        x_[0]                 = xmin;
        x_[N - 1]             = xmax;
        for (long i = 1; i < N - 1; ++i) {
            const double x2 = -1. + i * d2;  // x2 between -1 and 1;
            x_[i]           = midpoint + c1 * std::atan(c2 * std::tan(0.5 * M_PI * x2));
        }
    }

    min_   = std::min(xmin, xmax);
    max_   = std::max(xmin, xmax);
    start_ = xmin;
    end_   = xmax;
}

FocusSpacing::Spec FocusSpacing::spec() const {
    Spec spacing_specs;
    spacing_specs.set("type", static_type());
    spacing_specs.set("N", size());
    spacing_specs.set("start", start_);
    spacing_specs.set("end", end_);
    spacing_specs.set("focus_factor", focus_factor_);
    return spacing_specs;
}

namespace {
static SpacingBuilder<FocusSpacing> __builder(FocusSpacing::static_type());
}

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
