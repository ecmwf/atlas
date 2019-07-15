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

#include "DistributionImpl.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace grid {

namespace {
std::string distribution_type( int N, const Partitioner& p = Partitioner() ) {
    if ( N == 1 ) {
        return "serial";
    }
    if ( not p ) {
        return "custom";
    }
    return p.type();
}
}  // namespace

DistributionImpl::DistributionImpl( const Grid& grid ) :
    nb_partitions_( 1 ),
    part_( grid.size(), 0 ),
    nb_pts_( nb_partitions_, grid.size() ),
    max_pts_( grid.size() ),
    min_pts_( grid.size() ),
    type_( distribution_type( nb_partitions_ ) ) {}

DistributionImpl::DistributionImpl( const Grid& grid, const Partitioner& partitioner ) {
    part_.resize( grid.size() );
    partitioner.partition( grid, part_.data() );
    nb_partitions_ = partitioner.nb_partitions();
    nb_pts_.resize( nb_partitions_, 0 );
    for ( idx_t j = 0, size = static_cast<idx_t>( part_.size() ); j < size; ++j ) {
        ++nb_pts_[part_[j]];
    }
    max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
    min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
    type_    = distribution_type( nb_partitions_, partitioner );
}

DistributionImpl::DistributionImpl( idx_t npts, int part[], int part0 ) {
    part_.assign( part, part + npts );
    std::set<int> partset( part_.begin(), part_.end() );
    nb_partitions_ = static_cast<idx_t>( partset.size() );
    nb_pts_.resize( nb_partitions_, 0 );
    for ( idx_t j = 0, size = static_cast<idx_t>( part_.size() ); j < size; ++j ) {
        part_[j] -= part0;
        ++nb_pts_[part_[j]];
    }
    max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
    min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
    type_    = distribution_type( nb_partitions_ );
}

DistributionImpl::~DistributionImpl() {}

void DistributionImpl::print( std::ostream& s ) const {
    s << "Distribution( "
      << "type: " << type_ << ", nbPoints: " << part_.size() << ", nbPartitions: " << nb_pts_.size() << ", parts : [";
    for ( idx_t i = 0, size = static_cast<idx_t>( part_.size() ); i < size; i++ ) {
        if ( i != 0 ) {
            s << ',';
        }
        s << part_[i];
    }
    s << ']';
}


DistributionImpl* atlas__GridDistribution__new( idx_t npts, int part[], int part0 ) {
    return new DistributionImpl( npts, part, part0 );
}

void atlas__GridDistribution__delete( DistributionImpl* This ) {
    delete This;
}

}  // namespace grid
}  // namespace atlas
