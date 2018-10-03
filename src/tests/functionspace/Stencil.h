#include "eckit/memory/ScopedPtr.h"
#include "eckit/types/Types.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/MicroDeg.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {


class ComputeVertical {
    std::vector<idx_t> nvaux_;
    idx_t nlev_;
    idx_t nlevaux_;

public:
    ComputeVertical( const std::vector<double>& zcoord ) {
        nlev_          = zcoord.size() - 2;
        double dzcoord = zcoord.back() - zcoord.front();
        ASSERT( dzcoord > 0 );
        for ( idx_t jlev = 0; jlev <= nlev_; ++jlev ) {
            dzcoord = std::min( dzcoord, zcoord[jlev + 1] - zcoord[jlev] );
        }
        nlevaux_ = std::round( 2. / dzcoord + 0.5 ) + 1;
        nvaux_.resize( nlevaux_ + 1 );
        double dzaux = 1. / double( nlevaux_ );

        double zaux = 0.;
        for ( idx_t jlevaux = 0; jlevaux <= nlevaux_; ++jlevaux ) {
            for ( idx_t jlev = 0; jlev <= nlev_; ++jlev ) {
                if ( zaux <= zcoord[jlev + 1] ) {
                    nvaux_[jlevaux] = jlev;
                    break;
                }
            }
            zaux += dzaux;
        }
    }
    idx_t operator()( double z ) { return nvaux_[std::floor( z * nlevaux_ )]; }
};


class ComputeNorth {
    std::vector<double> y_;
    double dy_;
    constexpr double tol() const { return 0.5e-6; }

public:
    ComputeNorth( const grid::StructuredGrid& grid ) {
        y_.resize( grid.ny() );
        for ( idx_t j = 0; j < grid.ny(); ++j ) {
            y_[j] = grid.y( j ) + tol();
        }
        dy_ = std::abs( grid.y( 1 ) - grid.y( 0 ) );
    }
    idx_t operator()( double y ) const {
        idx_t j = std::floor( ( y_[0] - y ) / dy_ );
#ifndef NDEBUG
        ASSERT( j >= -1 );
#endif
        while ( y_[std::max( j, idx_t{0} )] > y ) {
            ++j;
        }
        if ( j >= 0 ) {
            do {
                --j;
            } while ( j >= y_.size() || y_[j] < y );
        }
        return j;
    }
};

class ComputeWest {
    std::vector<double> dx;
    std::vector<double> xref;
    constexpr double tol() const { return 0.5e-6; }
    idx_t ny;

public:
    ComputeWest( const grid::StructuredGrid& grid ) {
        if ( not grid::RegularGrid( grid ) &&
             std::abs( std::max( std::abs( grid.y().front() ), std::abs( grid.y().back() ) ) - 90. ) < tol() ) {
            throw eckit::NotImplemented( "ComputeWest not yet implemented for irregular grids with latitudes at pole",
                                         Here() );
        }
        dx.resize( grid.ny() );
        xref.resize( grid.ny() );
        for ( idx_t j = 0; j < grid.ny(); ++j ) {
            dx[j]   = std::abs( grid.x( 1, j ) - grid.x( 0, j ) );
            xref[j] = grid.x( 0, j ) - tol();
        }
        ny = grid.ny();
    }
    idx_t operator()( const double& x, idx_t j ) const {
        idx_t jj{j};
        if ( jj < 0 ) { jj = -j - 1; }
        if ( jj >= ny ) { jj = ny - 1 - ( jj - ny ); }
        idx_t i = std::floor( ( x - xref[jj] ) / dx[jj] );
        return i;
    }
};

class ComputeHorizontalStencil;

template <idx_t StencilWidth>
class HorizontalStencil {
    friend class ComputeHorizontalStencil;
    std::array<idx_t, StencilWidth> i_begin_;
    idx_t j_begin_;

public:
    idx_t i( idx_t offset ) const { return i_begin_[offset]; }
    idx_t j( idx_t offset ) const { return j_begin_ + offset; }
    constexpr idx_t width() const { return StencilWidth; }
};


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
    ComputeNorth compute_north_;
    ComputeWest compute_west_;
    idx_t stencil_width_;
    idx_t stencil_begin_;

public:
    ComputeHorizontalStencil( const grid::StructuredGrid& grid, idx_t stencil_width ) :
        compute_north_( grid ),
        compute_west_( grid ),
        stencil_width_( stencil_width ) {
        stencil_begin_ = stencil_width_ - idx_t( double( stencil_width_ ) / 2. + 1. );
    }
    template <typename Vector>
    void operator()( const double& x, const double& y, Vector& i, idx_t& j ) const {
        j = compute_north_( y ) - stencil_begin_;
        for ( idx_t jj = 0; jj < stencil_width_; ++jj ) {
            i[jj] = compute_west_( x, j + jj ) - stencil_begin_;
        }
    }
    HorizontalStencil<4> operator()( const double& x, const double& y ) const {
        HorizontalStencil<4> stencil;
        operator()( x, y, stencil );
        return stencil;
    }
    template <typename stencil_t>
    void operator()( const double& x, const double& y, stencil_t& stencil ) const {
        operator()( x, y, stencil.i_begin_, stencil.j_begin_ );
    }
};

}  // namespace test
}  // namespace atlas
