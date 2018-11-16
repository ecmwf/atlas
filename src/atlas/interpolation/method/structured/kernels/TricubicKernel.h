/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#include <cmath>
#include <limits>

#include "eckit/exception/Exceptions.h"
#include "eckit/linalg/Triplet.h"

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"

#include "BicubicKernel.h"

namespace atlas {
namespace interpolation {
namespace method {


class CubicVerticalInterpolation {
    ComputeVerticalStencil compute_vertical_stencil_;
    Vertical vertical_;
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }
    idx_t first_level_;
    idx_t last_level_;
    bool limiter_;

public:
    CubicVerticalInterpolation() = default;

    CubicVerticalInterpolation( const Vertical& vertical, const eckit::Configuration& config = util::NoConfig() ) :
        compute_vertical_stencil_( vertical, stencil_width() ),
        vertical_( vertical ),
        first_level_( vertical_.k_begin() ),
        last_level_( vertical_.k_end() - 1 ) {
        limiter_ = config.getBool( "limiter", false );
    }
    struct Weights {
        std::array<double, 4> weights_k;
    };
    using Stencil = VerticalStencil<4>;

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

        //        auto quadratic_interpolation = [z]( const double zvec[], double w[] ) {
        //            double d01 = zvec[0] - zvec[1];
        //            double d02 = zvec[0] - zvec[2];
        //            double d12 = zvec[1] - zvec[2];
        //            double dc0 = d01 * d02;
        //            double dc1 = -d01 * d12;
        //            double d0  = z - zvec[0];
        //            double d1  = z - zvec[1];
        //            double d2  = z - zvec[2];
        //            w[0]       = ( d1 * d2 ) / dc0;
        //            w[1]       = ( d0 * d2 ) / dc1;
        //            w[2]       = 1. - w[0] - w[1];
        //        };

        if ( stencil.k_interval() == -1 ) {
            // constant extrapolation
            //        lev0   lev1   lev2   lev3
            //      +  |------X------X------X
            //        w=1    w=0    w=0    w=0
            w[0] = 1.;
            w[1] = 0.;
            w[2] = 0.;
            w[3] = 0.;
            return;
        }
        else if ( stencil.k_interval() == 3 ) {
            // constant extrapolation
            //   lev(n-4)  lev(n-3)  lev(n-2)  lev(n-1)
            //      X---------X---------X---------|   +
            //     w=0      w=0       w=0       w=1
            w[0] = 0.;
            w[1] = 0.;
            w[2] = 0.;
            w[3] = 1.;
            return;
        }
        //        else if ( stencil.k_interval() == 0 ) {
        //            // quadratic interpolation
        //            // lev0   lev1   lev2   lev3
        //            //  |  +   |      |      |
        //            //                      w=0
        //            quadratic_interpolation( zvec.data(), w.data() );
        //            w[3] = 0.;
        //            return;
        //        }
        //        else if ( stencil.k_interval() == 2 ) {
        //            // quadratic interpolation
        //            // lev(n-4)  lev(n-3)  lev(n-2)  lev(n-1)
        //            //   |         |         |    +    |
        //            //  w=0
        //            quadratic_interpolation( zvec.data() + 1, w.data() + 1 );
        //            w[0] = 0.;
        //            return;
        //        }

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
            output += w[k] * input[stencil.k( k )];
        }


        if ( limiter_ ) {
            idx_t k = stencil.k_interval();
            idx_t k1, k2;
            if ( k < 0 ) { k1 = k2 = 0; }
            else if ( k > 2 ) {
                k1 = k2 = 3;
            }
            else {
                k1 = k;
                k2 = k + 1;
            }
            double f1     = input[stencil.k( k1 )];
            double f2     = input[stencil.k( k2 )];
            double maxval = std::max( f1, f2 );
            double minval = std::min( f1, f2 );
            output        = std::min( maxval, std::max( minval, output ) );
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

class TricubicKernel {
    using Triplet  = eckit::linalg::Triplet;
    using Triplets = std::vector<Triplet>;

public:
    TricubicKernel( const functionspace::StructuredColumns& fs ) {
        src_ = fs;
        ASSERT( src_ );
        ASSERT( src_.halo() >= 2 );
        horizontal_interpolation_ = BicubicKernel( src_ );
        vertical_interpolation_   = CubicVerticalInterpolation( fs.vertical() );
    }

private:
    functionspace::StructuredColumns src_;
    BicubicKernel horizontal_interpolation_;
    CubicVerticalInterpolation vertical_interpolation_;
    bool limiter_{false};

public:
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width() * stencil_width(); }
    static constexpr idx_t stencil_halo() {
        return static_cast<idx_t>( static_cast<double>( stencil_width() ) / 2. + 0.5 );
    }

public:
    using Stencil = Stencil3D<4>;
    struct Weights {
        std::array<std::array<double, 4>, 4> weights_i;
        std::array<double, 4> weights_j;
        std::array<double, 4> weights_k;
    };

public:
    struct WorkSpace {
        Stencil stencil;
        Weights weights;
    };

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
    typename std::enable_if<( array_t::RANK == 2 ), typename array_t::value_type>::type interpolate(
        const stencil_t& stencil, const weights_t& weights, const array_t& input ) const {
        using Value = typename array_t::value_type;

        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& wj = weights.weights_j;
        const auto& wk = weights.weights_k;

        Value output = 0.;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& wi = weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t n   = src_.index( stencil.i( i, j ), stencil.j( j ) );
                Value wij = wi[i] * wj[j];
                for ( idx_t k = 0; k < stencil_width(); ++k ) {
                    Value w = wij * wk[k];
                    output += w * input( n, stencil.k( k ) );
                }
                index[j][i] = n;
            }
        }

        if ( limiter_ ) { limit( output, index, stencil, input ); }
        return output;
    }

    template <typename array_t, typename stencil_t>
    typename std::enable_if<( array_t::RANK == 2 ), void>::type limit( typename array_t::value_type& output,
                                                                       const std::array<std::array<idx_t, 4>, 4>& index,
                                                                       const stencil_t& stencil,
                                                                       const array_t& input ) const {
        using Scalar = typename array_t::value_type;
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        idx_t k = stencil.k_interval();
        idx_t k1, k2;
        if ( k < 0 ) { k1 = k2 = 0; }
        else if ( k > 2 ) {
            k1 = k2 = 3;
        }
        else {
            k1 = k;
            k2 = k + 1;
        }

        Scalar maxval = std::numeric_limits<Scalar>::lowest();
        Scalar minval = std::numeric_limits<Scalar>::max();
        for ( idx_t j = 1; j < 3; ++j ) {
            for ( idx_t i = 1; i < 3; ++i ) {
                idx_t n = index[j][i];

                Scalar f1 = input( n, stencil.k( k1 ) );
                Scalar f2 = input( n, stencil.k( k2 ) );

                maxval = std::max( maxval, f1 );
                maxval = std::max( maxval, f2 );
                minval = std::min( minval, f1 );
                minval = std::min( minval, f2 );
            }
        }
        if ( output < minval ) { output = minval; }
        else if ( output > maxval ) {
            output = maxval;
        }
    }


    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<( OutputArray::RANK == 1 ), void>::type interpolate( const stencil_t& stencil,
                                                                                 const weights_t& weights,
                                                                                 const InputArray& input,
                                                                                 OutputArray& output, idx_t r ) const {
        output( r ) = interpolate( stencil, weights, input );
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<( OutputArray::RANK == 2 ), void>::type interpolate( const stencil_t& stencil,
                                                                                 const weights_t& weights,
                                                                                 const InputArray& input,
                                                                                 OutputArray& output, idx_t r,
                                                                                 idx_t k ) const {
        output( r, k ) = interpolate( stencil, weights, input );
    }

    template <typename array_t>
    typename array_t::value_type operator()( const double x, const double y, const double z,
                                             const array_t& input ) const {
        Stencil stencil;
        compute_stencil( x, y, z, stencil );
        Weights weights;
        compute_weights( x, y, z, stencil, weights );
        return interpolate( stencil, weights, input );
    }

    template <typename array_t>
    typename array_t::value_type interpolate( const PointLonLat& p, const double z, const array_t& input,
                                              WorkSpace& ws ) const {
        compute_stencil( p.lon(), p.lat(), z, ws.stencil );
        compute_weights( p.lon(), p.lat(), z, ws.stencil, ws.weights );
        return interpolate( ws.stencil, ws.weights, input );
    }

    // Thread private workspace
    Triplets compute_triplets( const idx_t row, const double x, const double y, const double z, WorkSpace& ws ) const {
        Triplets triplets;
        triplets.reserve( stencil_size() );
        insert_triplets( row, x, y, z, triplets, ws );
        return triplets;
    }

    Triplets reserve_triplets( size_t N ) {
        Triplets triplets;
        triplets.reserve( N * stencil_size() );
        return triplets;
    }

    Triplets allocate_triplets( size_t N ) { return Triplets( N * stencil_size() ); }

    void insert_triplets( const idx_t row, const PointXY& p, const double z, Triplets& triplets, WorkSpace& ws ) const {
        insert_triplets( row, p.x(), p.y(), z, triplets, ws );
    }

    void insert_triplets( const idx_t row, const double x, const double y, const double z, Triplets& triplets,
                          WorkSpace& ws ) const {
        compute_stencil( x, y, z, ws.stencil );
        compute_weights( x, y, z, ws.stencil, ws.weights );
        const auto& wj = ws.weights.weights_j;
        const auto& wk = ws.weights.weights_k;

        idx_t pos = row * stencil_size();
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& wi = ws.weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t n    = src_.index( ws.stencil.i( i, j ), ws.stencil.j( j ) );
                double wij = wi[i] * wj[j];
                for ( idx_t k = 0; k < stencil_width(); ++k ) {
                    idx_t col       = n * src_.vertical().size() + ws.stencil.k( k );
                    double w        = wij * wk[k];
                    triplets[pos++] = Triplet( row, col, w );
                }
            }
        }
    }
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
