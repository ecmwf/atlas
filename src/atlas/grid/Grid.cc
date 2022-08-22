/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/Grid.h"

#include <limits>
#include <string>
#include <vector>

#include "atlas/domain/Domain.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/grid/Spacing.h"
#include "atlas/projection/Projection.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"

namespace atlas {

Grid::IterateXY Grid::xy() const {
    return Grid::IterateXY(*get());
}

Grid::IterateLonLat Grid::lonlat() const {
    return Grid::IterateLonLat(*get());
}

Grid::Grid(const std::string& shortname, const Domain& domain):
    Handle([&] {
        Config config;
        if (domain) {
            config.set("domain", domain.spec());
        }
        return Grid::Implementation::create(shortname, config);
    }()) {}

Grid::Grid(const std::string& shortname, const Projection& projection, const Domain& domain):
    Handle([&] {
        Config config;
        if (projection) {
            config.set("projection", projection.spec());
        }
        if (domain) {
            config.set("domain", domain.spec());
        }
        return Grid::Implementation::create(shortname, config);
    }()) {}

Grid::Grid(const Grid& grid, const Grid::Domain& domain):
    Handle([&] {
        ATLAS_ASSERT(grid);
        return Grid::Implementation::create(*grid.get(), domain);
    }()) {}

Grid::Grid(const Config& p): Handle(Grid::Implementation::create(p)) {}

idx_t Grid::size() const {
    return get()->size();
}

size_t Grid::footprint() const {
    return get()->footprint();
}

const Grid::Projection& Grid::projection() const {
    return get()->projection();
}

const Grid::Domain& Grid::domain() const {
    return get()->domain();
}

RectangularLonLatDomain Grid::lonlatBoundingBox() const {
    return get()->lonlatBoundingBox();
}

std::string Grid::name() const {
    return get()->name();
}

std::string Grid::type() const {
    return get()->type();
}

std::string Grid::uid() const {
    return get()->uid();
}

void Grid::hash(eckit::Hash& h) const {
    get()->hash(h);
}

Grid::Spec Grid::spec() const {
    return get()->spec();
}

Grid::Config Grid::meshgenerator() const {
    return get()->meshgenerator();
}

Grid::Config Grid::partitioner() const {
    return get()->partitioner();
}

}  // namespace atlas
