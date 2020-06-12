/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "BandsPartitioner.h"

#include <cstdint>

#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/detail/distribution/BandsDistribution.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

size_t BandsPartitioner::blocksize( const atlas::grid::detail::partitioner::Partitioner::Grid& grid ) const {
    if ( blocksize_ == BLOCKSIZE_NX ) {
        if ( auto regular = RegularGrid( grid ) ) {
            return regular.nx();
        }
        return 1;
    }
    ATLAS_ASSERT( blocksize_ > 0 );
    return size_t( blocksize_ );
}

BandsPartitioner::BandsPartitioner() : Partitioner(), blocksize_{1} {}

BandsPartitioner::BandsPartitioner( int N, int blocksize ) : Partitioner( N ), blocksize_( blocksize ) {}

Distribution BandsPartitioner::partition( const Partitioner::Grid& grid ) const {
    if ( not distribution::BandsDistribution<int32_t>::detectOverflow( grid.size(), nb_partitions(),
                                                                       blocksize( ( grid ) ) ) ) {
        return Distribution{
            new distribution::BandsDistribution<int32_t>{grid, nb_partitions(), type(), blocksize( grid )}};
    }
    else {
        return Distribution{
            new distribution::BandsDistribution<int64_t>{grid, nb_partitions(), type(), blocksize( grid )}};
    }
}

void BandsPartitioner::partition( const Partitioner::Grid& grid, int part[] ) const {
    gidx_t gridsize = grid.size();
    if ( not distribution::BandsDistribution<int32_t>::detectOverflow( grid.size(), nb_partitions(),
                                                                       blocksize( ( grid ) ) ) ) {
        distribution::BandsDistribution<int32_t> distribution{grid, nb_partitions(), type(), blocksize( grid )};
        for ( gidx_t n = 0; n < gridsize; ++n ) {
            part[n] = distribution.function( n );
        }
    }
    else {
        distribution::BandsDistribution<int64_t> distribution{grid, nb_partitions(), type(), blocksize( grid )};
        for ( gidx_t n = 0; n < gridsize; ++n ) {
            part[n] = distribution.function( n );
        }
    }
}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas

namespace {
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::BandsPartitioner> __Bands(
    atlas::grid::detail::partitioner::BandsPartitioner::static_type() );
}
