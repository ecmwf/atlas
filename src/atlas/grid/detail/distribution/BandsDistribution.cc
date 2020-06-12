/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "BandsDistribution.h"

#include <algorithm>
#include <cstdint>

#include "atlas/grid/Grid.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {
namespace detail {
namespace distribution {

template <typename Int>
BandsDistribution<Int>::BandsDistribution( const atlas::Grid& grid, atlas::idx_t nb_partitions, const std::string& type,
                                           size_t blocksize ) :
    DistributionFunctionT<BandsDistribution<Int>>( grid ) {
    this->type_          = type;
    size_t gridsize      = grid.size();
    this->size_          = gridsize;
    this->nb_partitions_ = nb_partitions;
    nb_partitions_Int_   = nb_partitions;
    blocksize_           = blocksize;

    nb_blocks_ = gridsize / blocksize_;

    if ( gridsize % blocksize_ )
        nb_blocks_++;

    this->nb_pts_.reserve( nb_partitions_Int_ );

    for ( idx_t iproc = 0; iproc < nb_partitions; iproc++ ) {
        // Approximate values
        gidx_t imin = blocksize_ * ( ( ( iproc + 0 ) * nb_blocks_ ) / nb_partitions_Int_ );
        gidx_t imax = blocksize_ * ( ( ( iproc + 1 ) * nb_blocks_ ) / nb_partitions_Int_ );

        while ( imin > 0 ) {
            if ( function( imin - blocksize_ ) == iproc )
                imin -= blocksize_;
            else
                break;
        }

        while ( function( imin ) < iproc ) {
            imin += blocksize_;
        }

        while ( function( imax - 1 ) == iproc + 1 ) {
            imax -= blocksize_;
        }

        while ( imax + blocksize_ <= gridsize ) {
            if ( function( imax ) == iproc ) {
                imax += blocksize_;
            }
            else {
                break;
            }
        }

        imax = std::min( imax, (gidx_t)gridsize );
        this->nb_pts_.push_back( imax - imin );
    }

    this->max_pts_ = *std::max_element( this->nb_pts_.begin(), this->nb_pts_.end() );
    this->min_pts_ = *std::min_element( this->nb_pts_.begin(), this->nb_pts_.end() );

    ATLAS_ASSERT( detectOverflow( gridsize, nb_partitions_Int_, blocksize_ ) == false );
}

template <typename Int>
bool BandsDistribution<Int>::detectOverflow( size_t gridsize, size_t nb_partitions, size_t blocksize ) {
    int64_t size                 = gridsize;
    int64_t iblock               = ( size - 1 ) / int64_t( blocksize );
    int64_t intermediate_product = iblock * int64_t( nb_partitions );
    ATLAS_ASSERT( intermediate_product > 0, "Even 64 bits is insufficient to prevent overflow." );
    if ( intermediate_product > std::numeric_limits<Int>::max() ) {
        return true;
    }
    return false;
}


template class BandsDistribution<int32_t>;
template class BandsDistribution<int64_t>;

}  // namespace distribution
}  // namespace detail
}  // namespace grid
}  // namespace atlas
