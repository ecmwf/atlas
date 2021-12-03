/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/spacing/CustomSpacing.h"

#include <algorithm>

#include "eckit/config/Parametrisation.h"

#include "atlas/grid/detail/spacing/SpacingFactory.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {
namespace spacing {

CustomSpacing::CustomSpacing(long N, const double values[], const Interval& interval) {
    x_.assign(values, values + N);
    min_ = std::min(interval[0], interval[1]);
    max_ = std::max(interval[0], interval[1]);
}

CustomSpacing::CustomSpacing(const eckit::Parametrisation& params) {
    params.get("values", x_);

    size_t N;
    if (params.get("N", N)) {
        ATLAS_ASSERT(x_.size() == N);
    }
    N = x_.size();

    std::vector<double> interval;
    if (params.get("interval", interval)) {
        min_ = std::min(interval[0], interval[1]);
        max_ = std::max(interval[0], interval[1]);
    }
    else {
        min_ = x_.front();
        max_ = x_.front();
        for (size_t j = 1; j < N; ++j) {
            min_ = std::min(min_, x_[j]);
            max_ = std::max(max_, x_[j]);
        }
    }
}

CustomSpacing::Spec CustomSpacing::spec() const {
    Spec spacing_specs;
    spacing_specs.set("type", static_type());
    spacing_specs.set("values", x_);
    spacing_specs.set("interval", std::vector<double>{min(), max()});
    return spacing_specs;
}

namespace {
static SpacingBuilder<CustomSpacing> __builder(CustomSpacing::static_type());
}

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
