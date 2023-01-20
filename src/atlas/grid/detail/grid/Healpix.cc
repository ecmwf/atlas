/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Healpix.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>

#include "atlas/grid/Spacing.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/grid/detail/spacing/CustomSpacing.h"
#include "atlas/util/Constants.h"

namespace atlas {
namespace grid {
namespace {  // anonymous

//---------------------------------------------------------------------------------------------------------------------


static class HealpixGridBuilder : GridBuilder {
public:
    HealpixGridBuilder(): GridBuilder("healpix", {"^[Hh]([0-9]+)$"}, {"H<n>"}) {}

    void print(std::ostream& os) const override {
        os << std::left << std::setw(20) << "H<n>"
           << "Healpix grid with <n>x<n> points for each of 12 tiles";
    }

    const Grid::Implementation* create(const std::string& name, const Grid::Config& config) const override {
        int id;
        std::vector<std::string> matches;
        if (match(name, matches, id)) {
            int N = std::stoi(matches[0]);
            return new detail::grid::Healpix(N);
        }
        return nullptr;
    }

    const Grid::Implementation* create(const Grid::Config& config) const override {
        long N;
        config.get("N", N);
        return new detail::grid::Healpix(N);
    }

} healpix_builder_;

//---------------------------------------------------------------------------------------------------------------------

}  // namespace

namespace detail {
namespace grid {

Healpix::XSpace healpix_xspace(long N) {
    std::vector<Spacing> vp(4 * N - 1);

    // Polar caps
    for (int r = 1; r < N; r++) {
        double start      = 45./r;
        vp[r - 1]         = LinearSpacing(start, start + 360., 4 * r, false);
        vp[4 * N - r - 1] = vp[r - 1];
    }

    // Equatorial belt
    const double start = 45. / N;
    for (int r = N; r < 2 * N; r++) {
        double r_start    = start * (2. - (r - N + 1) % 2);
        vp[r - 1]         = LinearSpacing(r_start, r_start + 360., 4 * N, false);
        vp[4 * N - r - 1] = vp[r - 1];
    }

    // Equator
    double r_start = start * (1 - (N % 2 ? 1 : 0));
    vp[2 * N - 1]  = LinearSpacing(r_start, r_start + 360., 4 * N, false);

    return vp;
}

Healpix::YSpace healpix_yspace(long N) {
    constexpr double rad2deg = util::Constants::radiansToDegrees();
    std::vector<double> y(4 * N - 1);

    // Polar caps
    for (int r = 1; r < N; r++) {
        y[r - 1]         = 90. - rad2deg * std::acos(1. - r * r / (3. * N * N));
        y[4 * N - 1 - r] = -y[r - 1];
    }

    // Equatorial belt
    for (int r = N; r < 2 * N; r++) {
        y[r - 1]         = 90. - rad2deg * std::acos((4. * N - 2. * r) / (3. * N));
        y[4 * N - 1 - r] = -y[r - 1];
    }

    // Equator
    y[2 * N - 1] = 0.;

    return new spacing::CustomSpacing(y.size(), y.data());
}

Healpix::Healpix(long N):
    Structured("H" + std::to_string(N), healpix_xspace(N), healpix_yspace(N), Projection(), GlobalDomain()) {}


Healpix::Config Healpix::meshgenerator() const {
    return util::Config("type", "healpix");
}
Healpix::Config Healpix::partitioner() const {
    return util::Config("type", "equal_regions");
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
