/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "DistributionImpl.h"

#include <algorithm>

#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/distribution/DistributionArray.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {

DistributionImpl* atlas__GridDistribution__new( idx_t size, int part[], int part0 ) {
    return new detail::distribution::DistributionArray( 0, size, part, part0 );
}

void atlas__GridDistribution__delete( DistributionImpl* This ) {
    delete This;
}

void atlas__GridDistribution__nb_pts( DistributionImpl* This, idx_t nb_pts[] ) {
    const auto& nb_pts_ = This->nb_pts();
    std::copy( nb_pts_.begin(), nb_pts_.end(), &nb_pts[0] );
}

idx_t atlas__atlas__GridDistribution__nb_partitions( DistributionImpl* This ) {
    return This->nb_partitions();
}

int atlas__GridDistribution__partition_int32( DistributionImpl* dist, int i ) {
    return dist->partition( i );
}

int atlas__GridDistribution__partition_int64( DistributionImpl* dist, long i ) {
    return dist->partition( i );
}

DistributionImpl* atlas__GridDistribution__new__Grid_Config( const detail::grid::Grid* grid,
                                                             const eckit::Parametrisation* config ) {
    ATLAS_ASSERT( grid != nullptr, "grid is an invalid pointer" );
    ATLAS_ASSERT( config != nullptr, "config is an invalid pointer" );
    DistributionImpl* distribution;
    {
        Distribution d{Grid{grid}, *config};
        distribution = d.get();
        distribution->attach();
    }
    distribution->detach();
    return distribution;
}

DistributionImpl* atlas__GridDistribution__new__Grid_Partitioner(
    const detail::grid::Grid* grid, const detail::partitioner::Partitioner* partitioner ) {
    ATLAS_ASSERT( grid != nullptr, "grid is an invalid pointer" );
    ATLAS_ASSERT( partitioner != nullptr, "partitioner is an invalid pointer" );
    DistributionImpl* distribution;
    {
        Distribution d{Grid{grid}, Partitioner{partitioner}};
        distribution = d.get();
        distribution->attach();
    }
    distribution->detach();
    return distribution;
}


}  // namespace grid
}  // namespace atlas
