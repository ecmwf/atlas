/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <ostream>

#include "atlas/domain/detail/DomainFactory.h"
#include "atlas/domain/detail/ZonalBandDomain.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace domain {

namespace {

static bool _is_global(double ymin, double ymax) {
    const double eps = 1.e-12;
    return std::abs((ymax - ymin) - 180.) < eps;
}

static std::array<double, 2> get_interval_y(const eckit::Parametrisation& params) {
    constexpr double invalid = std::numeric_limits<double>::max();
    auto is_valid            = [](double y) {
        //  deepcode ignore FloatingPointEquals: We want exact comparison
        return (y != invalid);
    };
    double ymax = invalid;
    double ymin = invalid;

    if (params.get("south", ymin) || params.get("north", ymax)) {
        ymax = is_valid(ymax) ? ymax : 90.;
        ymin = is_valid(ymin) ? ymin : -90.;
    }
    else {
        if (!params.get("ymin", ymin)) {
            throw_Exception("ymin missing in Params", Here());
        }
        if (!params.get("ymax", ymax)) {
            throw_Exception("ymax missing in Params", Here());
        }
    }

    return {ymin, ymax};
}

static double get_west(const eckit::Parametrisation& params) {
    double west = 0.;
    params.get("west", west);
    return west;
}

/*
constexpr std::array<double, 2> interval_x() {
    return {0., 360.};
}
*/
}  // namespace

constexpr char ZonalBandDomain::units_[];

bool ZonalBandDomain::is_global(const Interval& y) {
    const double eps = 1.e-12;
    return std::abs(std::abs(y[1] - y[0]) - 180.) < eps;
}

ZonalBandDomain::ZonalBandDomain(const eckit::Parametrisation& params):
    ZonalBandDomain(get_interval_y(params), get_west(params)) {}

ZonalBandDomain::ZonalBandDomain(const Interval& interval_y): ZonalBandDomain(interval_y, /*west*/ 0.) {}

ZonalBandDomain::ZonalBandDomain(const Interval& interval_y, const double west):
    RectangularLonLatDomain({west, west + 360.}, interval_y) {
    global_   = _is_global(ymin(), ymax());
    ymin_tol_ = ymin() - 1.e-6;
    ymax_tol_ = ymax() + 1.e-6;
}

bool ZonalBandDomain::contains(double /*x*/, double y) const {
    return contains_y(y);
}

ZonalBandDomain::Spec ZonalBandDomain::spec() const {
    Spec domain_spec;
    domain_spec.set("type", type());
    domain_spec.set("ymin", ymin());
    domain_spec.set("ymax", ymax());
    if (xmin() != 0.) {
        domain_spec.set("west", xmin());
    }
    return domain_spec;
}

void ZonalBandDomain::hash(eckit::Hash& h) const {
    spec().hash(h);
}

void ZonalBandDomain::print(std::ostream& os) const {
    os << "ZonalBandDomain["
       << "ymin=" << ymin() << ",ymax=" << ymax() << ",west=" << xmin() << "]";
}

bool ZonalBandDomain::containsNorthPole() const {
    return ymax_tol_ >= 90.;
}

bool ZonalBandDomain::containsSouthPole() const {
    return ymin_tol_ <= -90.;
}

namespace {
static DomainBuilder<ZonalBandDomain> register_builder(ZonalBandDomain::static_type());
}

}  // namespace domain
}  // namespace atlas
