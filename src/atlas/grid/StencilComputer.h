/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <vector>
#include "atlas/grid/Grid.h"
#include "atlas/grid/Vertical.h"
#include "atlas/library/config.h"

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------
/// @class ComputeVertical
/// @brief Compute lower vertical level index for given coordinate
/// zcoord:
/// @verbatim
///   0----1----2----3--...--(n-1)----(n)----(n+1)
///   --->|<---|<---|<--...-|<--------------------
/// @endverbatim
/// If coordinate falls on vertical level (+- epsilon), that level is returned
/// If coordinate falls in range [0,1) or [n,n+1],
/// the index is snapped to 1 and (n-1) respectively. This allows reliably that
/// the returned index can be used for stencil operations.
///
/// IFS full levels don't have a level at the boundaries (0.,1.)
/// It is the "half" levels that contain (0.,1.). For reasons of boundary conditions
/// however, the full levels also have 0. prepended and 1. appended.
///
/// Example IFS full levels for regular distribution dz ( level 0 and n+1 are added for boundary conditions )
///  0      :  0.0
///  jlev   :  jlev*dz - 0.5*dz
///  nlev   :  nlev*dz - 0.5*dz
///  nlev+1 :  1.0

class ComputeLower {
    std::vector<double> z_;
    std::vector<idx_t> nvaux_;
    idx_t nlev_;
    idx_t nlevaux_;
    double rlevaux_;

public:
    ComputeLower( const Vertical& z );

    idx_t operator()( double z ) const {
        idx_t idx = static_cast<idx_t>( std::floor( z * rlevaux_ ) );
#ifndef NDEBUG
        ASSERT( idx < static_cast<idx_t>( nvaux_.size() ) && idx >= 0 );
#endif
        idx = nvaux_[idx];
        if ( idx < nlev_ - 1 && z > z_[idx + 1] ) { ++idx; }
        return idx;
    }
};

//-----------------------------------------------------------------------------

class ComputeNorth {
    std::vector<double> y_;
    double dy_;
    static constexpr double tol() { return 0.5e-6; }
    idx_t halo_;
    idx_t ny_;

public:
    ComputeNorth( const grid::StructuredGrid& grid, idx_t halo );

    idx_t operator()( double y ) const {
        idx_t j = static_cast<idx_t>( std::floor( ( y_[halo_ + 0] - y ) / dy_ ) );
        while ( y_[halo_ + j] > y ) {
            ++j;
        }
        do {
            --j;
        } while ( y_[halo_ + j] < y );

        return j;
    }
};

//-----------------------------------------------------------------------------

class ComputeWest {
    std::vector<double> dx;
    std::vector<double> xref;
    idx_t halo_;  // halo in north-south direction
    idx_t ny_;
    static constexpr double tol() { return 0.5e-6; }

public:
    ComputeWest( const grid::StructuredGrid& grid, idx_t halo = 0 );

    idx_t operator()( const double& x, idx_t j ) const {
        idx_t jj = halo_ + j;
        idx_t i  = static_cast<idx_t>( std::floor( ( x - xref[jj] ) / dx[jj] ) );
        return i;
    }
};  // namespace test


//-----------------------------------------------------------------------------

// @class ComputeHorizontalStencil
// @brief Compute stencil in horizontal direction (i,j)
//
// Given a stencil width, the stencil for a given P{x,y} is:
//
//        i[0]     i[1]     i[2]    i[3]
//         x        x        x         x       j + 0
//          x       x       x        x         j + 1
//                     P
//          x       x       x        x         j + 2
//         x        x        x         x       j + 3
//
//   In case the x-component of P is aligned with any
//   stencil, gridpoint, the stencil will assume the grid point
//   is on the point P's left side:
//
//        i[0]     i[1]     i[2]    i[3]
//         x        x        x         x       j + 0
//          x       x       x        x         j + 1
//                  P
//          x       x       x        x         j + 2
//         x        x        x         x       j + 3

class ComputeHorizontalStencil {
    idx_t halo_;
    ComputeNorth compute_north_;
    ComputeWest compute_west_;
    idx_t stencil_width_;
    idx_t stencil_begin_;

public:
    ComputeHorizontalStencil( const grid::StructuredGrid& grid, idx_t stencil_width );

    template <typename stencil_t>
    void operator()( const double& x, const double& y, stencil_t& stencil ) const {
        stencil.j_begin_ = compute_north_( y ) - stencil_begin_;
        for ( idx_t jj = 0; jj < stencil_width_; ++jj ) {
            stencil.i_begin_[jj] = compute_west_( x, stencil.j_begin_ + jj ) - stencil_begin_;
        }
    }
};


//-----------------------------------------------------------------------------

class ComputeVerticalStencil {
    ComputeLower compute_lower_;
    idx_t stencil_width_;
    idx_t stencil_begin_;

public:
    ComputeVerticalStencil( const Vertical& vertical, idx_t stencil_width );

    template <typename stencil_t>
    void operator()( const double& z, stencil_t& stencil ) const {
        stencil.k_begin_ = compute_lower_( z ) - stencil_begin_;
    }
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas