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
#include "eckit/linalg/SparseMatrix.h"

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
    static constexpr double tol() { return 0.5e-6; }
    static constexpr double halo() { return 5; }
    idx_t ny_;

public:
    ComputeNorth( const grid::StructuredGrid& grid ) {
        ny_ = grid.ny();
        y_.resize( ny_ + 2 * halo() );
        ASSERT( halo() < ny_ );
        idx_t north_pole_included = 90. - std::abs( grid.y().front() ) < tol();
        idx_t south_pole_included = 90. - std::abs( grid.y().back() ) < tol();

        for ( idx_t j = -halo(); j < 0; ++j ) {
            idx_t jj       = -j - 1 + north_pole_included;
            y_[halo() + j] = 180. - grid.y( jj ) + tol();
        }
        for ( idx_t j = 0; j < ny_; ++j ) {
            y_[halo() + j] = grid.y( j ) + tol();
        }
        for ( idx_t j = ny_; j < ny_ + halo(); ++j ) {
            idx_t jj       = 2 * ny_ - j - 1 - south_pole_included;
            y_[halo() + j] = -180. - grid.y( jj ) + tol();
        }
        dy_ = std::abs( grid.y( 1 ) - grid.y( 0 ) );
    }

    idx_t operator()( double y ) const {
        idx_t j = std::floor( ( y_[halo() + 0] - y ) / dy_ );
        while ( y_[halo() + j] > y ) {
            ++j;
        }
        do {
            --j;
        } while ( y_[halo() + j] < y );

        return j;
    }
};

class ComputeWest {
    std::vector<double> dx;
    std::vector<double> xref;
    static constexpr double tol() { return 0.5e-6; }
    static constexpr double halo() { return 5; }
    idx_t ny_;

public:
    ComputeWest( const grid::StructuredGrid& grid ) {
        idx_t north_pole_included = 90. - std::abs( grid.y().front() ) < tol();
        idx_t south_pole_included = 90. - std::abs( grid.y().back() ) < tol();
        ny_                       = grid.ny();
        dx.resize( ny_ + 2 * halo() );
        xref.resize( ny_ + 2 * halo() );
        for ( idx_t j = -halo(); j < 0; ++j ) {
            idx_t jj         = -j - 1 + north_pole_included;
            dx[halo() + j]   = grid.x( 1, jj ) - grid.x( 0, jj );
            xref[halo() + j] = grid.x( 0, jj ) - tol();
        }
        for ( idx_t j = 0; j < ny_; ++j ) {
            dx[halo() + j]   = std::abs( grid.x( 1, j ) - grid.x( 0, j ) );
            xref[halo() + j] = grid.x( 0, j ) - tol();
        }
        for ( idx_t j = ny_; j < ny_ + halo(); ++j ) {
            idx_t jj         = 2 * ny_ - j - 1 - south_pole_included;
            dx[halo() + j]   = std::abs( grid.x( 1, jj ) - grid.x( 0, jj ) );
            xref[halo() + j] = grid.x( 0, jj ) - tol();
        }
    }
    idx_t operator()( const double& x, idx_t j ) const {
        idx_t jj = halo() + j;
        idx_t i  = std::floor( ( x - xref[jj] ) / dx[jj] );
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

class CubicStructuredInterpolation {
    ComputeHorizontalStencil compute_horizontal_stencil_;
    functionspace::StructuredColumns fs_;
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width()*stencil_width(); }
    bool limiter_{true};

public:
    CubicStructuredInterpolation( const functionspace::StructuredColumns& fs ) :
        fs_( fs ),
        compute_horizontal_stencil_( fs.grid(), stencil_width() ) {}
    template <typename weights_t>
    void compute_weights( const double x, const double y, weights_t& weights ) const {
        HorizontalStencil<stencil_width()> stencil;
        compute_horizontal_stencil_( x, y, stencil );
        compute_weights( x, y, stencil, weights );
    }

    struct Weights {
        std::array<std::array<double, 4>, 4> weights_i_foreach_j;
        std::array<double, 4> weights_j;
    };

    template <typename stencil_t, typename weights_t>
    void compute_weights( const double x, const double y, const stencil_t& stencil, weights_t& weights ) const {
        PointXY P1, P2;
        std::array<double, 4> yvec;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            auto& weights_i = weights.weights_i_foreach_j[j];
            fs_.compute_xy( stencil.i( j ) + 1, stencil.j( j ), P1 );
            fs_.compute_xy( stencil.i( j ) + 2, stencil.j( j ), P2 );
            double alpha     = ( P2.x() - x ) / ( P2.x() - P1.x() );
            double alpha_sqr = alpha * alpha;
            weights_i[0]     = -alpha * ( 1. - alpha_sqr ) / 6.;
            weights_i[1]     = 0.5 * alpha * ( 1. + alpha ) * ( 2. - alpha );
            weights_i[2]     = 0.5 * ( 1. - alpha_sqr ) * ( 2. - alpha );
            weights_i[3]     = 1. - weights_i[0] - weights_i[1] - weights_i[2];
            yvec[j]          = P1.y();
        }
        double dl12 = yvec[0] - yvec[1];
        double dl13 = yvec[0] - yvec[2];
        double dl14 = yvec[0] - yvec[3];
        double dl23 = yvec[1] - yvec[2];
        double dl24 = yvec[1] - yvec[3];
        double dl34 = yvec[2] - yvec[3];
        double dcl1 = dl12 * dl13 * dl14;
        double dcl2 = -dl12 * dl23 * dl24;
        double dcl3 = dl13 * dl23 * dl34;

        double dl1 = y - yvec[0];
        double dl2 = y - yvec[1];
        double dl3 = y - yvec[2];
        double dl4 = y - yvec[3];

        auto& weights_j = weights.weights_j;
        weights_j[0]    = ( dl2 * dl3 * dl4 ) / dcl1;
        weights_j[1]    = ( dl1 * dl3 * dl4 ) / dcl2;
        weights_j[2]    = ( dl1 * dl2 * dl4 ) / dcl3;
        weights_j[3]    = 1. - weights_j[0] - weights_j[1] - weights_j[2];
    }

    template <typename stencil_t, typename weights_t, typename array_t>
    void interpolate( const stencil_t& stencil, const weights_t& weights, const array_t& input, double& output ) {
        std::array<double, 4> output_j;
        std::fill( std::begin( output_j ), std::end( output_j ), 0 );
        double maxval = std::numeric_limits<double>::lowest();
        double minval = std::numeric_limits<double>::max();
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& weights_i = weights.weights_i_foreach_j[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t n    = fs_.index( stencil.i( j ) + i, stencil.j( j ) );
                double inp = input( n );
                output_j[j] += weights_i[i] * inp;
                maxval = std::max( maxval, inp );
                minval = std::min( minval, inp );
            }
        }
        output          = 0.;
        const auto& weights_j = weights.weights_j;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            output += weights_j[j] * output_j[j];
        }
        if ( limiter_ ) { output = std::min( maxval, std::max( minval, output ) ); }
    }

    template <typename array_t>
    double operator()( const double x, const double y, const array_t& input ) {
        HorizontalStencil<stencil_width()> stencil;
        compute_horizontal_stencil_( x, y, stencil );
        Weights weights;
        compute_weights( x, y, stencil, weights );
        double output;
        interpolate( stencil, weights, input, output );
        return output;
    }

    struct WorkSpace {
        HorizontalStencil<4> stencil;
        Weights weights;
    };

    using Triplet  = eckit::linalg::Triplet;
    using Triplets = std::vector<Triplet>;

    // Thread private workspace
    Triplets compute_triplets( const idx_t row, const double x, const double y, WorkSpace& ws ) {
        Triplets triplets; triplets.reserve( stencil_size() );
        insert_triplets( row, x, y, triplets, ws );
        return triplets;
    }

    Triplets reserve_triplets( size_t N ) {
        Triplets triplets;
        triplets.reserve( N * stencil_size() );
        return triplets;
    }

    void insert_triplets( const idx_t row, const double x, const double y, Triplets& triplets, WorkSpace& ws ) {
        compute_horizontal_stencil_( x, y, ws.stencil );
        compute_weights( x, y, ws.stencil, ws.weights );
        const auto& wj = ws.weights.weights_j;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& wi = ws.weights.weights_i_foreach_j[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t col = fs_.index( ws.stencil.i( j ) + i, ws.stencil.j( j ) );
                double w = wi[i]*wj[j];
                triplets.emplace_back( row, col, w );
            }
        }
    }
};

}  // namespace test
}  // namespace atlas
