/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>

#include "Distribution.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/distribution/DistributionImpl.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace grid {

Distribution::Distribution( const Grid& grid ) : Handle( new Implementation( grid ) ) {}
Distribution::Distribution( const Grid& grid, const Config & config ) : Handle( new Implementation( grid, config ) ) {}

Distribution::Distribution( const Grid& grid, const Partitioner& partitioner ) :
    Handle( new Implementation( grid, partitioner ) ) {}

Distribution::Distribution( int nb_partitions, idx_t npts, int part[], int part0 ) :
    Handle( new Implementation( nb_partitions, npts, part, part0 ) ) {}

Distribution::Distribution( int nb_partitions, partition_t&& part ) :
    Handle( new Implementation( nb_partitions, std::move( part ) ) ) {}

Distribution::~Distribution() = default;

const Distribution::partition_t& Distribution::partition() const {
    return get()->partition();
}

idx_t Distribution::nb_partitions() const {
    return get()->nb_partitions();
}

const int* Distribution::data() const {
    return get()->data();
}

const std::vector<idx_t>& Distribution::nb_pts() const {
    return get()->nb_pts();
}

idx_t Distribution::max_pts() const {
    return get()->max_pts();
}

idx_t Distribution::min_pts() const {
    return get()->min_pts();
}

const std::string& Distribution::type() const {
    return get()->type();
}

std::ostream& operator<<( std::ostream& os, const Distribution& distribution ) {
    distribution.get()->print( os );
    return os;
}

Distribution::operator const partition_t&() const {
    return *get();
}

}  // namespace grid
}  // namespace atlas
