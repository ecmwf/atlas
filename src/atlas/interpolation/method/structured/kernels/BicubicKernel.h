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

namespace atlas {
namespace interpolation {
namespace method {

class BicubicKernel {
    using Triplet  = eckit::linalg::Triplet;
    using Triplets = std::vector<Triplet>;

public:
    BicubicKernel() = default;

    BicubicKernel( const functionspace::StructuredColumns& fs ) {
        src_ = fs;
        ASSERT( src_ );
        ASSERT( src_.halo() >= 2 );
        compute_horizontal_stencil_ = ComputeHorizontalStencil( src_.grid(), stencil_width() );
    }

private:
    functionspace::StructuredColumns src_;
    ComputeHorizontalStencil compute_horizontal_stencil_;
    bool limiter_{false};

public:
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }
    static constexpr idx_t stencil_halo() {
        return static_cast<idx_t>( static_cast<double>( stencil_width() ) / 2. + 0.5 );
    }

public:
    using Stencil = HorizontalStencil<4>;
    struct Weights {
        std::array<std::array<double, 4>, 4> weights_i;
        std::array<double, 4> weights_j;
    };

public:
    struct WorkSpace {
        Stencil stencil;
        Weights weights;
    };

    template <typename stencil_t>
    void compute_stencil( const double x, const double y, stencil_t& stencil ) const {
        compute_horizontal_stencil_( x, y, stencil );
    }

    template <typename weights_t>
    void compute_weights( const double x, const double y, weights_t& weights ) const {
        Stencil stencil;
        compute_stencil( x, y, stencil );
        compute_weights( x, y, stencil, weights );
    }


    template <typename stencil_t, typename weights_t>
    void compute_weights( const double x, const double y, const stencil_t& stencil, weights_t& weights ) const {
        PointXY P1, P2;
        std::array<double, 4> yvec;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            auto& weights_i = weights.weights_i[j];
            src_.compute_xy( stencil.i( 1, j ), stencil.j( j ), P1 );
            src_.compute_xy( stencil.i( 2, j ), stencil.j( j ), P2 );
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
    typename array_t::value_type interpolate( const stencil_t& stencil, const weights_t& weights,
                                              const array_t& input ) const {
        using Value = typename array_t::value_type;

        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& weights_j = weights.weights_j;
        Value output          = 0.;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& weights_i = weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t n = src_.index( stencil.i( i, j ), stencil.j( j ) );
                Value w = weights_i[i] * weights_j[j];
                output += w * input[n];
                index[j][i] = n;
            }
        }

        if ( limiter_ ) { limit( output, index, input ); }
        return output;
    }

    template <typename array_t>
    void limit( typename array_t::value_type& output, const std::array<std::array<idx_t, 4>, 4>& index,
                const array_t& input ) const {
        using Scalar = typename array_t::value_type;
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        Scalar maxval = std::numeric_limits<Scalar>::lowest();
        Scalar minval = std::numeric_limits<Scalar>::max();
        for ( idx_t j = 1; j < 3; ++j ) {
            for ( idx_t i = 1; i < 3; ++i ) {
                idx_t n    = index[j][i];
                Scalar val = input[n];
                maxval     = std::max( maxval, val );
                minval     = std::min( minval, val );
            }
        }
        if ( output < minval ) { output = minval; }
        else if ( output > maxval ) {
            output = maxval;
        }
    }


    template <typename stencil_t, typename weights_t, typename Value, int Rank>
    typename std::enable_if<( Rank == 1 ), void>::type interpolate( const stencil_t& stencil, const weights_t& weights,
                                                                    const array::ArrayView<Value, Rank>& input,
                                                                    array::ArrayView<Value, Rank>& output,
                                                                    idx_t r ) const {
        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& weights_j = weights.weights_j;
        output( r )           = 0.;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& weights_i = weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t n = src_.index( stencil.i( i, j ), stencil.j( j ) );
                Value w = static_cast<Value>( weights_i[i] * weights_j[j] );
                output( r ) += w * input[n];
                index[j][i] = n;
            }
        }

        if ( limiter_ ) { limit( index, input, output, r ); }
    }

    template <typename Value, int Rank>
    typename std::enable_if<( Rank == 1 ), void>::type limit( const std::array<std::array<idx_t, 4>, 4>& index,
                                                              const array::ArrayView<Value, Rank>& input,
                                                              array::ArrayView<Value, Rank>& output, idx_t r ) const {
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        Value maxval = std::numeric_limits<Value>::lowest();
        Value minval = std::numeric_limits<Value>::max();
        for ( idx_t j = 1; j < 3; ++j ) {
            for ( idx_t i = 1; i < 3; ++i ) {
                idx_t n   = index[j][i];
                Value val = input[n];
                maxval    = std::max( maxval, val );
                minval    = std::min( minval, val );
            }
        }
        if ( output( r ) < minval ) { output( r ) = minval; }
        else if ( output( r ) > maxval ) {
            output( r ) = maxval;
        }
    }


    template <typename stencil_t, typename weights_t, typename Value, int Rank>
    typename std::enable_if<( Rank == 2 ), void>::type interpolate( const stencil_t& stencil, const weights_t& weights,
                                                                    const array::ArrayView<Value, Rank>& input,
                                                                    array::ArrayView<Value, Rank>& output,
                                                                    idx_t r ) const {
        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& weights_j = weights.weights_j;
        const idx_t Nk        = output.shape( 1 );
        for ( idx_t k = 0; k < Nk; ++k ) {
            output( r, k ) = 0.;
        }
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& weights_i = weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t n = src_.index( stencil.i( i, j ), stencil.j( j ) );
                Value w = static_cast<Value>( weights_i[i] * weights_j[j] );
                for ( idx_t k = 0; k < Nk; ++k ) {
                    output( r, k ) += w * input( n, k );
                }
                index[j][i] = n;
            }
        }

        if ( limiter_ ) { limit( index, input, output, r ); }
    }

    template <typename Value, int Rank>
    typename std::enable_if<( Rank == 2 ), void>::type limit( const std::array<std::array<idx_t, 4>, 4>& index,
                                                              const array::ArrayView<Value, Rank>& input,
                                                              array::ArrayView<Value, Rank>& output, idx_t r ) const {
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        for ( idx_t k = 0; k < output.shape( 1 ); ++k ) {
            Value maxval = std::numeric_limits<Value>::lowest();
            Value minval = std::numeric_limits<Value>::max();
            for ( idx_t j = 1; j < 3; ++j ) {
                for ( idx_t i = 1; i < 3; ++i ) {
                    idx_t n   = index[j][i];
                    Value val = input( n, k );
                    maxval    = std::max( maxval, val );
                    minval    = std::min( minval, val );
                }
            }
            if ( output( r, k ) < minval ) { output( r, k ) = minval; }
            else if ( output( r, k ) > maxval ) {
                output( r, k ) = maxval;
            }
        }
    }


    template <typename array_t>
    typename array_t::value_type operator()( const double x, const double y, const array_t& input ) const {
        Stencil stencil;
        compute_stencil( x, y, stencil );
        Weights weights;
        compute_weights( x, y, stencil, weights );
        return interpolate( stencil, weights, input );
    }

    template <typename array_t>
    typename array_t::value_type interpolate( const PointLonLat& p, const array_t& input, WorkSpace& ws ) const {
        compute_stencil( p.lon(), p.lat(), ws.stencil );
        compute_weights( p.lon(), p.lat(), ws.stencil, ws.weights );
        return interpolate( ws.stencil, ws.weights, input );
    }

    // Thread private workspace
    Triplets compute_triplets( const idx_t row, const double x, const double y, WorkSpace& ws ) const {
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

    Triplets allocate_triplets( size_t N ) { return Triplets( N * stencil_size() ); }

    void insert_triplets( const idx_t row, const PointXY& p, Triplets& triplets, WorkSpace& ws ) const {
        insert_triplets( row, p.x(), p.y(), triplets, ws );
    }

    void insert_triplets( const idx_t row, const double x, const double y, Triplets& triplets, WorkSpace& ws ) const {
        compute_horizontal_stencil_( x, y, ws.stencil );
        compute_weights( x, y, ws.stencil, ws.weights );
        const auto& wj = ws.weights.weights_j;

        idx_t pos = row * stencil_size();
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& wi = ws.weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t col       = src_.index( ws.stencil.i( i, j ), ws.stencil.j( j ) );
                double w        = wi[i] * wj[j];
                triplets[pos++] = Triplet( row, col, w );
            }
        }
    }
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas