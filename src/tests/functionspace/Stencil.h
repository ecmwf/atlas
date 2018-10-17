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

#include "eckit/types/Types.h"
#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

class Vertical {
    idx_t k_begin_;
    idx_t k_end_;
    idx_t size_;
    bool boundaries_;
    std::vector<double> z_;

public:
    idx_t k_begin() const { return k_begin_; }
    idx_t k_end() const { return k_end_; }
    bool boundaries() const { return boundaries_; }
    idx_t size() const { return size_; }

    template <typename Int>
    double operator()( const Int k ) const {
        return z_[k];
    }

    template <typename Int>
    double operator[]( const Int k ) const {
        return z_[k];
    }

    template <typename Vector>
    Vertical( idx_t levels, const Vector& z, const util::Config& config = util::NoConfig() ) {
        size_       = levels;
        boundaries_ = config.getBool( "boundaries", false );
        k_begin_    = 0;
        k_end_      = size_;
        if ( boundaries_ ) {
            size_ += 2;
            ++k_begin_;
            k_end_ = size_ - 1;
        }
        ASSERT( size_ == z.size() );
        z_.resize( size_ );
        for ( idx_t k = 0; k < size_; ++k ) {
            z_[k] = z[k];
        }
    }

    const std::vector<double>& zcoord() const { return z_; }
};

// @class ComputeVertical
// @brief Helper to compute vertical level below
//
//
//  IFS full levels for regular distribution ( level 0 and nlev+1 are added for boundary conditions )
//  0      :  0.0
//  jlev   :  jlev*dz - 0.5*dz
//  nlev   :  nlev*dz - 0.5*dz
//  nlev+1 :  1.0

class ComputeVertical {
    std::vector<double> zcoord_;
    std::vector<idx_t> nvaux_;
    idx_t nlev_;
    idx_t nlevaux_;
    double rlevaux_;

public:
    template <typename Vector>
    ComputeVertical( const Vector& zcoord ) {
        nlev_ = zcoord.size() - 2;
        zcoord_.resize( nlev_ + 2 );
        double dzcoord       = std::numeric_limits<double>::max();
        constexpr double tol = 1.e-12;
        ASSERT( dzcoord > 0 );
        for ( idx_t jlev = 0; jlev < nlev_; ++jlev ) {
            dzcoord       = std::min( dzcoord, zcoord[jlev + 1] - zcoord[jlev] );
            zcoord_[jlev] = zcoord[jlev] - tol;
        }
        zcoord_[nlev_]     = zcoord_[nlev_] - tol;
        zcoord_[nlev_ + 1] = zcoord_[nlev_ + 1] - tol;

        nlevaux_ = static_cast<idx_t>( std::round( 2. / dzcoord + 0.5 ) + 1 );
        rlevaux_ = double( nlevaux_ );
        nvaux_.resize( nlevaux_ + 1 );
        double dzaux = ( zcoord[nlev_ + 1] - zcoord[0] ) / rlevaux_;

        idx_t iref = 1;
        for ( idx_t jlevaux = 0; jlevaux <= nlevaux_; ++jlevaux ) {
            if ( jlevaux * dzaux >= zcoord[iref + 1] && iref < nlev_ - 1 ) { ++iref; }
            nvaux_[jlevaux] = iref;
        }
    }

    idx_t operator()( double z ) const {
        idx_t idx = static_cast<idx_t>( std::floor( z * rlevaux_ ) );
#ifndef NDEBUG
        ASSERT( idx < static_cast<idx_t>( nvaux_.size() ) && idx >= 0 );
#endif
        // ATLAS_DEBUG_VAR( idx );
        idx = nvaux_[idx];
        // ATLAS_DEBUG_VAR( z );
        // ATLAS_DEBUG_VAR( zcoord_[idx] );
        // ATLAS_DEBUG_VAR( idx );
        if ( idx < nlev_ - 1 && z > zcoord_[idx + 1] ) {
            ++idx;
            // ATLAS_DEBUG();
        }
        return idx;
        //return nvaux_[idx];
    }
};


class ComputeNorth {
    std::vector<double> y_;
    double dy_;
    static constexpr double tol() { return 0.5e-6; }
    static constexpr idx_t halo() { return 5; }
    idx_t ny_;

public:
    ComputeNorth( const grid::StructuredGrid& grid ) {
        ASSERT( grid );
        if ( not grid.domain().global() ) {
            throw eckit::NotImplemented( "Only implemented for global grids", Here() );
        }

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
        idx_t j = static_cast<idx_t>( std::floor( ( y_[halo() + 0] - y ) / dy_ ) );
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
    static constexpr idx_t halo() { return 5; }
    idx_t ny_;

public:
    ComputeWest( const grid::StructuredGrid& grid ) {
        ASSERT( grid );
        if ( not grid.domain().global() ) {
            throw eckit::NotImplemented( "Only implemented for global grids", Here() );
        }
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
        idx_t i  = static_cast<idx_t>( std::floor( ( x - xref[jj] ) / dx[jj] ) );
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

template <idx_t StencilWidth>
class VerticalStencil {
    friend class ComputeVerticalStencil;
    idx_t k_begin_;

public:
    idx_t k( idx_t offset ) const { return k_begin_ + offset; }
    constexpr idx_t width() const { return StencilWidth; }
};

class ComputeVerticalStencil {
    ComputeVertical compute_vertical_;
    idx_t stencil_width_;
    idx_t stencil_begin_;

public:
    template <typename Vector>
    ComputeVerticalStencil( const Vector& zcoord, idx_t stencil_width ) :
        compute_vertical_( zcoord ),
        stencil_width_( stencil_width ) {
        stencil_begin_ = stencil_width_ - idx_t( double( stencil_width_ ) / 2. + 1. );
    }

    void operator()( const double& z, idx_t& k ) const { k = compute_vertical_( z ) - stencil_begin_; }
    VerticalStencil<4> operator()( const double& z ) const {
        VerticalStencil<4> stencil;
        operator()( z, stencil );
        return stencil;
    }
    template <typename stencil_t>
    void operator()( const double& z, stencil_t& stencil ) const {
        operator()( z, stencil.k_begin_ );
    }
};

class CubicVerticalInterpolation {
    ComputeVerticalStencil compute_vertical_stencil_;
    Vertical vertical_;
    bool boundaries_;
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }
    bool limiter_{false};
    idx_t first_level_;
    idx_t last_level_;

public:
    CubicVerticalInterpolation( const Vertical& vertical ) :
        compute_vertical_stencil_( vertical, stencil_width() ),
        vertical_( vertical ),
        boundaries_( vertical.boundaries() ),
        first_level_( vertical_.k_begin() ),
        last_level_( vertical_.k_end() - 1 ) {}
    struct Weights {
        std::array<double, 4> weights_k;
    };

    template <typename stencil_t>
    void compute_stencil( const double z, stencil_t& stencil ) const {
        compute_vertical_stencil_( z, stencil );
    }

    template <typename stencil_t, typename weights_t>
    void compute_weights( const double z, const stencil_t& stencil, weights_t& weights ) const {
        auto& w = weights.weights_k;

        std::array<double, 4> zvec;
        for ( idx_t k = 0; k < 4; ++k ) {
            zvec[k] = vertical_( stencil.k( k ) );
        }

        if ( boundaries_ ) {
            auto quadratic_interpolation = [z]( const double zvec[], double w[] ) {
                double d01 = zvec[0] - zvec[1];
                double d02 = zvec[0] - zvec[2];
                double d12 = zvec[1] - zvec[2];
                double dc0 = d01 * d02;
                double dc1 = -d01 * d12;
                double d0  = z - zvec[0];
                double d1  = z - zvec[1];
                double d2  = z - zvec[2];
                w[0]       = ( d1 * d2 ) / dc0;
                w[1]       = ( d0 * d2 ) / dc1;
                w[2]       = 1. - w[0] - w[1];
            };

            if ( z < vertical_( first_level_ ) or z > vertical_( last_level_ ) ) {
                // linear extrapolation
                // lev0   lev1   lev2   lev3            lev(n-2)  lev(n-1) lev(n)  lev(n+1)
                //  X   +  |      |      X      or         X        |        |  +    X
                //  w=0                 w=0               w=0                       w=0
                w[3] = 0.;
                w[2] = ( z - zvec[1] ) / ( zvec[2] - zvec[1] );
                w[1] = 1. - w[2];
                w[0] = 0.;
                return;
            }
            else if ( z < vertical_( first_level_ + 1 ) ) {
                // quadratic interpolation
                // lev0   lev1   lev2   lev3
                //  X      |  +   |      |
                //  w=0
                quadratic_interpolation( zvec.data() + 1, w.data() + 1 );
                w[0] = 0.;
                return;
            }
            else if ( z > vertical_( last_level_ - 1 ) ) {
                // quadratic interpolation
                // lev(n-2)  lev(n-1) lev(n)  lev(n+1)
                //   |        |   +    |       X
                //                            w=0
                quadratic_interpolation( zvec.data(), w.data() );
                w[3] = 0.;
                return;
            }
        }

        // cubic interpolation
        // lev(k+0)   lev(k+1)   lev(k+2)   lev(k+3)
        //    |          |     x    |          |
        double d01 = zvec[0] - zvec[1];
        double d02 = zvec[0] - zvec[2];
        double d03 = zvec[0] - zvec[3];
        double d12 = zvec[1] - zvec[2];
        double d13 = zvec[1] - zvec[3];
        double d23 = zvec[2] - zvec[3];
        double dc0 = d01 * d02 * d03;
        double dc1 = -d01 * d12 * d13;
        double dc2 = d02 * d12 * d23;

        double d0 = z - zvec[0];
        double d1 = z - zvec[1];
        double d2 = z - zvec[2];
        double d3 = z - zvec[3];

        w[0] = ( d1 * d2 * d3 ) / dc0;
        w[1] = ( d0 * d2 * d3 ) / dc1;
        w[2] = ( d0 * d1 * d3 ) / dc2;
        w[3] = 1. - w[0] - w[1] - w[2];
    }

    template <typename stencil_t, typename weights_t, typename array_t>
    void interpolate( const stencil_t& stencil, const weights_t& weights, const array_t& input, double& output ) const {
        output        = 0.;
        const auto& w = weights.weights_k;
        for ( idx_t k = 0; k < stencil_width(); ++k ) {
            output += w[k] * input( stencil.k( k ) );
        }
        if ( limiter_ ) {
            double f1     = input( stencil.k( 1 ) );
            double f2     = input( stencil.k( 2 ) );
            double maxval = std::max( f1, f2 );
            double minval = std::min( f1, f2 );
            ;
            output = std::min( maxval, std::max( minval, output ) );
        }
    }

    template <typename array_t>
    double operator()( const double z, const array_t& input ) const {
        VerticalStencil<stencil_width()> stencil;
        compute_vertical_stencil_( z, stencil );
        Weights weights;
        compute_weights( z, stencil, weights );
        double output;
        interpolate( stencil, weights, input, output );
        return output;
    }
};

class CubicStructuredInterpolation {
    functionspace::StructuredColumns fs_;
    ComputeHorizontalStencil compute_horizontal_stencil_;
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }
    bool limiter_{false};

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
        std::array<std::array<double, 4>, 4> weights_i;
        std::array<double, 4> weights_j;
    };

    template <typename stencil_t>
    void compute_stencil( const double x, const double y, stencil_t& stencil ) const {
        compute_horizontal_stencil_( x, y, stencil );
    }

    template <typename stencil_t, typename weights_t>
    void compute_weights( const double x, const double y, const stencil_t& stencil, weights_t& weights ) const {
        PointXY P1, P2;
        std::array<double, 4> yvec;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            auto& weights_i = weights.weights_i[j];
            fs_.compute_xy( stencil.i( j ) + 1, stencil.j( j ), P1 );
            fs_.compute_xy( stencil.i( j ) + 2, stencil.j( j ), P2 );
            double alpha               = ( P2.x() - x ) / ( P2.x() - P1.x() );
            double alpha_sqr           = alpha * alpha;
            double two_minus_alpha     = 2. - alpha;
            double one_minus_alpha_sqr = 1. - alpha_sqr;
            weights_i[0]               = -alpha * one_minus_alpha_sqr / 6.;
            weights_i[1]               = 0.5 * alpha * ( 1. + alpha ) * two_minus_alpha;
            weights_i[2]               = 0.5 * one_minus_alpha_sqr * two_minus_alpha;
            weights_i[3]               = 1. - weights_i[0] - weights_i[1] - weights_i[2];
            yvec[j]                    = P1.y();
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
        std::array<idx_t, stencil_size()> index;
        std::fill( std::begin( output_j ), std::end( output_j ), 0 );
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& weights_i = weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t n = fs_.index( stencil.i( j ) + i, stencil.j( j ) );
                output_j[j] += weights_i[i] * input( n );
                index[stencil_width() * j + i] = n;
            }
        }
        output                = 0.;
        const auto& weights_j = weights.weights_j;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            output += weights_j[j] * output_j[j];
        }
        if ( limiter_ ) {
            // Limit output to max/min of values in stencil marked by '*'
            //         x        x        x         x
            //              x     *-----*     x
            //                   /   P  |
            //          x       *------ *        x
            //        x        x        x         x
            double maxval = std::numeric_limits<double>::lowest();
            double minval = std::numeric_limits<double>::max();
            for ( idx_t j = 1; j < 3; ++j ) {
                for ( idx_t i = 1; i < 3; ++i ) {
                    idx_t n    = index[stencil_width() * j + i];
                    double val = input( n );
                    maxval     = std::max( maxval, val );
                    minval     = std::min( minval, val );
                }
            }
            output = std::min( maxval, std::max( minval, output ) );
        }
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
        Triplets triplets;
        triplets.reserve( stencil_size() );
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
            const auto& wi = ws.weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t col = fs_.index( ws.stencil.i( j ) + i, ws.stencil.j( j ) );
                double w  = wi[i] * wj[j];
                triplets.emplace_back( row, col, w );
            }
        }
    }
};


template <idx_t StencilWidth>
class Stencil3D {
    friend class ComputeHorizontalStencil;
    friend class ComputeVerticalStencil;
    std::array<idx_t, StencilWidth> i_begin_;
    idx_t j_begin_;
    idx_t k_begin_;

public:
    idx_t i( idx_t offset ) const { return i_begin_[offset]; }
    idx_t j( idx_t offset ) const { return j_begin_ + offset; }
    idx_t k( idx_t offset ) const { return k_begin_ + offset; }
    constexpr idx_t width() const { return StencilWidth; }
};

class Cubic3DInterpolation {
    functionspace::StructuredColumns fs_;
    CubicStructuredInterpolation horizontal_interpolation_;
    CubicVerticalInterpolation vertical_interpolation_;
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }

public:
    struct Weights {
        std::array<std::array<double, 4>, 4> weights_i;
        std::array<double, 4> weights_j;
        std::array<double, 4> weights_k;
    };

    using Stencil = Stencil3D<4>;

    Cubic3DInterpolation( const functionspace::StructuredColumns& fs, const Vertical& vertical ) :
        fs_( fs ),
        horizontal_interpolation_( fs.grid() ),
        vertical_interpolation_( vertical ) {}

    template <typename stencil_t>
    void compute_stencil( const double x, const double y, const double z, stencil_t& stencil ) const {
        horizontal_interpolation_.compute_stencil( x, y, stencil );
        vertical_interpolation_.compute_stencil( z, stencil );
    }

    template <typename weights_t>
    void compute_weights( const double x, const double y, const double z, weights_t& weights ) const {
        Stencil stencil;
        compute_stencil( x, y, z, stencil );
        compute_weights( x, y, z, stencil, weights );
    }

    template <typename stencil_t, typename weights_t>
    void compute_weights( const double x, const double y, const double z, const stencil_t& stencil,
                          weights_t& weights ) const {
        horizontal_interpolation_.compute_weights( x, y, stencil, weights );
        vertical_interpolation_.compute_weights( z, stencil, weights );
    }

    template <typename stencil_t, typename weights_t, typename array_t>
    void interpolate( const stencil_t& stencil, const weights_t& weights, const array_t& input, double& output ) {
        output         = 0.;
        const auto& wj = weights.weights_j;
        const auto& wk = weights.weights_k;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& wi = weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t n = fs_.index( stencil.i( j ) + i, stencil.j( j ) );
                for ( idx_t k = 0; k < stencil_width(); ++k ) {
                    output += wi[i] * wj[j] * wk[k] * input( n, stencil.k( k ) );
                }
            }
        }
    }

    template <typename array_t>
    double operator()( const double x, const double y, const double z, const array_t& input ) {
        Stencil stencil;
        Weights weights;
        compute_stencil( x, y, z, stencil );
        compute_weights( x, y, z, stencil, weights );
        double output;
        interpolate( stencil, weights, input, output );
        return output;
    }
};

}  // namespace test
}  // namespace atlas
