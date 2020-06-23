/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "DistributionArray.h"

#include <algorithm>
#include <ostream>
#include <vector>

#include "eckit/types/Types.h"
#include "eckit/utils/Hash.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {
namespace detail {
namespace distribution {

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

DistributionArray::DistributionArray( const Grid& grid, const eckit::Parametrisation& config ) :
    DistributionArray( grid, Partitioner{config} ) {}


DistributionArray::DistributionArray( const Grid& grid, const Partitioner& partitioner ) {
    part_.resize( grid.size() );
    partitioner.partition( grid, part_.data() );
    nb_partitions_ = partitioner.nb_partitions();

    size_t size     = part_.size();
    int num_threads = atlas_omp_get_max_threads();

    std::vector<std::vector<int> > nb_pts_per_thread( num_threads, std::vector<int>( nb_partitions_ ) );
    atlas_omp_parallel {
        int thread   = atlas_omp_get_thread_num();
        auto& nb_pts = nb_pts_per_thread[thread];
        atlas_omp_for( size_t j = 0; j < size; ++j ) {
            int p = part_[j];
            ++nb_pts[p];
        }
    }

    nb_pts_.resize( nb_partitions_, 0 );
    for ( int thread = 0; thread < num_threads; ++thread ) {
        for ( int p = 0; p < nb_partitions_; ++p ) {
            nb_pts_[p] += nb_pts_per_thread[thread][p];
        }
    }

    max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
    min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
    type_    = distribution_type( nb_partitions_, partitioner );
}

DistributionArray::DistributionArray( int nb_partitions, idx_t npts, int part[], int part0 ) {
    part_.assign( part, part + npts );
    if ( nb_partitions == 0 ) {
        std::set<int> partset( part_.begin(), part_.end() );
        nb_partitions_ = static_cast<idx_t>( partset.size() );
    }
    else {
        nb_partitions_ = nb_partitions;
    }
    nb_pts_.resize( nb_partitions_, 0 );
    for ( idx_t j = 0, size = static_cast<idx_t>( part_.size() ); j < size; ++j ) {
        part_[j] -= part0;
        ++nb_pts_[part_[j]];
    }
    max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
    min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
    type_    = distribution_type( nb_partitions_ );
}

DistributionArray::DistributionArray( int nb_partitions, partition_t&& part ) :
    nb_partitions_( nb_partitions ), part_( std::move( part ) ), nb_pts_( nb_partitions_, 0 ) {
    size_t size     = part_.size();
    int num_threads = atlas_omp_get_max_threads();
    std::vector<std::vector<int> > nb_pts_per_thread( num_threads, std::vector<int>( nb_partitions_ ) );
    atlas_omp_parallel {
        int thread   = atlas_omp_get_thread_num();
        auto& nb_pts = nb_pts_per_thread[thread];
        atlas_omp_for( size_t j = 0; j < size; ++j ) {
            int p = part_[j];
            ++nb_pts[p];
        }
    }
    for ( int thread = 0; thread < num_threads; ++thread ) {
        for ( int p = 0; p < nb_partitions_; ++p ) {
            nb_pts_[p] += nb_pts_per_thread[thread][p];
        }
    }

    max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
    min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
    type_    = distribution_type( nb_partitions_ );
}

DistributionArray::~DistributionArray() = default;

void DistributionArray::print( std::ostream& s ) const {
    auto print_partition = [&]( std::ostream& s ) {
        eckit::output_list<int> list_printer( s );
        for ( size_t i = 0, size = part_.size(); i < size; i++ ) {
            list_printer.push_back( part_[i] );
        }
    };
    s << "Distribution( "
      << "type: " << type_ << ", nb_points: " << size() << ", nb_partitions: " << nb_pts_.size() << ", parts : ";
    print_partition( s );
}

void DistributionArray::hash( eckit::Hash& hash ) const {
    for ( size_t i = 0; i < part_.size(); i++ ) {
        hash.add( part_[i] );
    }
}


}  // namespace distribution
}  // namespace detail
}  // namespace grid
}  // namespace atlas
