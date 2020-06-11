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

#include "atlas/grid/Grid.h"

namespace atlas {
namespace grid {
namespace detail {
namespace distribution {

BandsDistribution::BandsDistribution( const atlas::Grid& grid, atlas::idx_t nb_partitions, const std::string& type,
                                      size_t blocksize ) :
    DistributionFunctionT<BandsDistribution>( grid ) {
    type_           = type;
    size_t gridsize = grid.size();
    size_           = gridsize;
    nb_partitions_  = nb_partitions;
    blocksize_      = blocksize;

    nb_blocks_ = gridsize / blocksize_;

    if ( gridsize % blocksize_ )
        nb_blocks_++;

    nb_pts_.reserve( nb_partitions_ );

    for ( idx_t iproc = 0; iproc < nb_partitions_; iproc++ ) {
        // Approximate values
        gidx_t imin = blocksize_ * ( ( ( iproc + 0 ) * nb_blocks_ ) / nb_partitions_ );
        gidx_t imax = blocksize_ * ( ( ( iproc + 1 ) * nb_blocks_ ) / nb_partitions_ );

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
        nb_pts_.push_back( imax - imin );
    }

    max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
    min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
}

}  // namespace distribution
}  // namespace detail
}  // namespace grid
}  // namespace atlas
