/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/config/Parametrisation.h"

#include "atlas/grid/detail/spacing/GaussianSpacing.h"
#include "atlas/grid/detail/spacing/SpacingFactory.h"
#include "atlas/grid/detail/spacing/gaussian/Latitudes.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {
namespace spacing {

GaussianSpacing::GaussianSpacing(long N) {
    // perform checks
    ATLAS_ASSERT(N % 2 == 0);

    // initialize latitudes during setup, to avoid repeating it.
    x_.resize(N);
    gaussian::gaussian_latitudes_npole_spole(N / 2, x_.data());

    min_ = -90.;
    max_ = 90.;
}

GaussianSpacing::GaussianSpacing(const eckit::Parametrisation& params) {
    // retrieve N from params
    long N;
    if (!params.get("N", N)) {
        throw_Exception("N missing in Params", Here());
    }

    // perform checks
    ATLAS_ASSERT(N % 2 == 0);

    // initialize latitudes during setup, to avoid repeating it.
    x_.resize(N);
    gaussian::gaussian_latitudes_npole_spole(N / 2, x_.data());

    // Not yet implemented: specify different bounds or direction (e.g from south
    // to north pole)
    double start = 90.;
    double end   = -90.;
    params.get("start", start);
    params.get("end", end);

    std::vector<double> interval;
    if (params.get("interval", interval)) {
        start = interval[0];
        end   = interval[1];
    }
    if (start != 90. && end != -90.) {
        ATLAS_NOTIMPLEMENTED;
    }

    min_ = std::min(start, end);
    max_ = std::max(start, end);
}

GaussianSpacing::Spec GaussianSpacing::spec() const {
    Spec spacing_specs;
    spacing_specs.set("type", static_type());
    spacing_specs.set("N", size());
    return spacing_specs;
}

namespace {
static SpacingBuilder<GaussianSpacing> __builder(GaussianSpacing::static_type());
}

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
