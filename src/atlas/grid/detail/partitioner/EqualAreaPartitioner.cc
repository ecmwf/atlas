/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/partitioner/EqualAreaPartitioner.h"

#include "atlas/grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/util/Constants.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

EqualAreaPartitioner::EqualAreaPartitioner():
    Partitioner(), partitioner_() {
}

EqualAreaPartitioner::EqualAreaPartitioner(int N):
    EqualAreaPartitioner(N, util::NoConfig()) {
}

EqualAreaPartitioner::EqualAreaPartitioner(int N, const eckit::Parametrisation& config):
    Partitioner(N,config), partitioner_(N,config) {
}

EqualAreaPartitioner::EqualAreaPartitioner(const eckit::Parametrisation& config):
    Partitioner(config), partitioner_(config) {
}

void EqualAreaPartitioner::partition(const Grid& grid, int part[]) const {

    if( partitioner_.coordinates_ == EqualRegionsPartitioner::Coordinates::XY && StructuredGrid(grid) ) {
        StructuredGrid g(grid);
        size_t n = 0;
        for (idx_t j=0; j<g.ny(); ++j) {
            const double lat = g.y(j) * util::Constants::degreesToRadians();
            int b = partitioner_.band(lat);
            int p = 0;
            for (int k = 0; k < b; ++k) {
                p += partitioner_.nb_regions(k);
            }
            idx_t nx = g.nx(j);
            for (idx_t i=0; i<nx; ++i) {
                const double lon = g.x(i,j) * util::Constants::degreesToRadians();
                part[n++] = p + partitioner_.sector(b, lon);
            }
        }
    }
    else {
        size_t j{0};
        for (PointLonLat p : grid.lonlat()) {
            p.lon() *= util::Constants::degreesToRadians();
            p.lat() *= util::Constants::degreesToRadians();
            part[j++] = partitioner_.partition(p.lon(), p.lat());
        }
    }

}


}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas

namespace {
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::EqualAreaPartitioner>
    __EqualArea("equal_area");
}
