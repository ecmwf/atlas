/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/StencilComputer.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace grid {

ComputeLower::ComputeLower( const Vertical& z ) : nlev_{z.size()} {
    z_.resize( nlev_ );
    double dz            = std::numeric_limits<double>::max();
    constexpr double tol = 1.e-12;
    ATLAS_ASSERT( dz > 0 );
    for ( idx_t jlev = 0; jlev < nlev_; ++jlev ) {
        if ( jlev + 1 < nlev_ ) {
            dz = std::min( dz, z[jlev + 1] - z[jlev] );
        }
        z_[jlev] = z[jlev] - tol;
    }

    nlevaux_ = static_cast<idx_t>( std::round( 2. * ( z.max() - z.min() ) / dz + 0.5 ) + 1 );
    rlevaux_ = double( nlevaux_ );
    nvaux_.resize( nlevaux_ + 1 );
    double dzaux = ( z.max() - z.min() ) / rlevaux_;

    idx_t iref = 0;
    for ( idx_t jlevaux = 0; jlevaux <= nlevaux_; ++jlevaux ) {
        if ( iref + 1 < nlev_ && jlevaux * dzaux >= z[iref + 1] ) {
            ++iref;
        }
        nvaux_[jlevaux] = iref;
    }
}

ComputeNorth::ComputeNorth( const StructuredGrid& grid, idx_t halo ) {
    ATLAS_ASSERT( grid );
    if ( not grid.domain().global() ) {
        throw_NotImplemented( "Only implemented for global grids", Here() );
    }
    halo_ = halo;
    ny_   = grid.ny();
    y_.resize( ny_ + 2 * halo_ );
    ATLAS_ASSERT( halo_ < ny_ );
    idx_t north_pole_included = 90. - std::abs( grid.y().front() ) < tol();
    idx_t south_pole_included = 90. - std::abs( grid.y().back() ) < tol();

    for ( idx_t j = -halo_; j < 0; ++j ) {
        idx_t jj      = -j - 1 + north_pole_included;
        y_[halo_ + j] = 180. - grid.y( jj ) + tol();
    }
    for ( idx_t j = 0; j < ny_; ++j ) {
        y_[halo_ + j] = grid.y( j ) + tol();
    }
    for ( idx_t j = ny_; j < ny_ + halo_; ++j ) {
        idx_t jj      = 2 * ny_ - j - 1 - south_pole_included;
        y_[halo_ + j] = -180. - grid.y( jj ) + tol();
    }
    dy_ = std::abs( grid.y( 1 ) - grid.y( 0 ) );
}

ComputeWest::ComputeWest( const StructuredGrid& grid, idx_t halo ) {
    ATLAS_ASSERT( grid );
    if ( not grid.domain().global() ) {
        throw_NotImplemented( "Only implemented for global grids", Here() );
    }
    halo_                     = halo;
    idx_t north_pole_included = 90. - std::abs( grid.y().front() ) < tol();
    idx_t south_pole_included = 90. - std::abs( grid.y().back() ) < tol();
    ny_                       = grid.ny();
    dx.resize( ny_ + 2 * halo_ );
    xref.resize( ny_ + 2 * halo_ );
    for ( idx_t j = -halo_; j < 0; ++j ) {
        idx_t jj        = -j - 1 + north_pole_included;
        dx[halo_ + j]   = grid.x( 1, jj ) - grid.x( 0, jj );
        xref[halo_ + j] = grid.x( 0, jj ) - tol();
    }
    for ( idx_t j = 0; j < ny_; ++j ) {
        dx[halo_ + j]   = std::abs( grid.x( 1, j ) - grid.x( 0, j ) );
        xref[halo_ + j] = grid.x( 0, j ) - tol();
    }
    for ( idx_t j = ny_; j < ny_ + halo_; ++j ) {
        idx_t jj        = 2 * ny_ - j - 1 - south_pole_included;
        dx[halo_ + j]   = std::abs( grid.x( 1, jj ) - grid.x( 0, jj ) );
        xref[halo_ + j] = grid.x( 0, jj ) - tol();
    }
}

ComputeHorizontalStencil::ComputeHorizontalStencil( const StructuredGrid& grid, idx_t stencil_width ) :
    halo_( ( stencil_width + 1 ) / 2 ),
    compute_north_( grid, halo_ ),
    compute_west_( grid, halo_ ),
    stencil_width_( stencil_width ) {
    stencil_begin_ = stencil_width_ - idx_t( double( stencil_width_ ) / 2. + 1. );
}

ComputeVerticalStencil::ComputeVerticalStencil( const Vertical& vertical, idx_t stencil_width ) :
    compute_lower_( vertical ), stencil_width_( stencil_width ) {
    stencil_begin_ = stencil_width_ - idx_t( double( stencil_width_ ) / 2. + 1. );
    clip_begin_    = 0;
    clip_end_      = vertical.size();
    vertical_min_  = vertical[clip_begin_];
    vertical_max_  = vertical[clip_end_ - 1];
}


//---------------------------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
