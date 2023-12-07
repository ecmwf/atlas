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

#include "atlas/grid/Grid.h"
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

void EqualAreaPartitioner::partition(const Grid& g, int part[]) const {
    size_t j{0};
    for (PointLonLat p : g.lonlat()) {
        p.lon() *= util::Constants::degreesToRadians();
        p.lat() *= util::Constants::degreesToRadians();
        part[j++] = partitioner_.partition(p.lon(), p.lat());
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
