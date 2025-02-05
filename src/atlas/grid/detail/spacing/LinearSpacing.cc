/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/spacing/LinearSpacing.h"

#include <cmath>
#include <ostream>


#include "eckit/config/Parametrisation.h"

#include "atlas/grid/detail/spacing/SpacingFactory.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {
namespace spacing {

LinearSpacing::Params::Params(const eckit::Parametrisation& params) {
    endpoint = true;
    params.get("endpoint", endpoint);
    std::vector<double> interval;
    if (params.get("N", N)) {
        // Only one remaining combinations possible:
        if (params.get("start", start) && params.get("end", end)) {
            // OK
        }
        else if (params.get("interval", interval)) {
            start = interval[0];
            end   = interval[1];
        }
        else if (params.get("start", start) && params.get("length", length)) {
            end = start + length;
        }
        else {
            throw_Exception("Invalid combination of parameters", Here());
        }
    }
    else {
        throw_Exception("Invalid combination of parameters", Here());
    }
    length = end - start;

    if (endpoint && N > 1) {
        step = length / double(N - 1);
    }
    else {
        step = length / double(N);
    }
}

LinearSpacing::Params::Params(double _start, double _end, long _N, bool _endpoint):
    start(_start), end(_end), N(_N), endpoint(_endpoint) {
    length = end - start;
    if (endpoint && N > 1) {
        step = length / double(N - 1);
    }
    else {
        step = length / double(N);
    }
}

LinearSpacing::LinearSpacing(const eckit::Parametrisation& params) {
    Params p(params);
    setup(p.start, p.end, p.N, p.endpoint);
}

LinearSpacing::LinearSpacing(double start, double end, long N, bool endpoint) {
    setup(start, end, N, endpoint);
}

LinearSpacing::LinearSpacing(const Interval& interval, long N, bool endpoint) {
    setup(interval[0], interval[1], N, endpoint);
}

void LinearSpacing::setup(double start, double end, long N, bool endpoint) {
    x_.resize(N);

    double step;
    volatile double _N = N;  // volatile keyword prevents agressive optimization by Cray compiler
    //  deepcode ignore FloatingPointEquals: Expect possible bit-identical start and end
    if (start == end) {
        step = 0.;
    }
    else if (endpoint && N > 1) {
        step = (end - start) / (_N - 1);
    }
    else {
        step = (end - start) / _N;
    }

    for (long i = 0; i < N; ++i) {
        x_[i] = start + i * step;
    }

    min_ = std::min(start, end);
    max_ = std::max(start, end);

    start_    = start;
    end_      = end;
    N_        = N;
    endpoint_ = endpoint;

    // For exact comparisons:
    if (N > 1) {
        x_.front() = start;
    }
    if (N > 2 && endpoint) {
        x_.back() = end;
    }
}

double LinearSpacing::step() const {
    if (size() > 1) {
        return x_[1] - x_[0];
    }
    else {
        return 0.;
    }
}

bool LinearSpacing::endpoint() const {
    return std::abs(x_.back() - max_) < 1.e-12;
}

LinearSpacing::Spec LinearSpacing::spec() const {
    Spec spacing_specs;
    spacing_specs.set("type", static_type());
    spacing_specs.set("start", start_);
    spacing_specs.set("end", end_);
    spacing_specs.set("N", N_);
    spacing_specs.set("endpoint", endpoint_);
    return spacing_specs;
}

namespace {
static SpacingBuilder<LinearSpacing> __builder(LinearSpacing::static_type());
}

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
