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

#include "eckit/utils/Hash.h"

#include "atlas/domain/detail/DomainFactory.h"
#include "atlas/domain/detail/GlobalDomain.h"

namespace atlas {
namespace domain {

namespace {
constexpr std::array<double, 2> yrange() {
    return {-90., 90.};
}

static double get_west(const eckit::Parametrisation& params) {
    double west = 0.;
    params.get("west", west);
    return west;
}

}  // namespace

GlobalDomain::GlobalDomain(const double west): ZonalBandDomain(yrange(), west) {}

GlobalDomain::GlobalDomain(): ZonalBandDomain(yrange()) {}

GlobalDomain::GlobalDomain(const eckit::Parametrisation& params): GlobalDomain(get_west(params)) {}

GlobalDomain::Spec GlobalDomain::spec() const {
    Spec domain_spec;
    domain_spec.set("type", type());
    if (xmin() != 0.) {
        domain_spec.set("west", xmin());
    }
    return domain_spec;
}

void GlobalDomain::hash(eckit::Hash& h) const {
    h.add(type());
}

void GlobalDomain::print(std::ostream& os) const {
    os << "GlobalDomain[west=" << xmin() << "]";
}

namespace {
static DomainBuilder<GlobalDomain> register_builder(GlobalDomain::static_type());
}
}  // namespace domain
}  // namespace atlas
