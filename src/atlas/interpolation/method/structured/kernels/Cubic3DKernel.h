/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction. and Interpolation
 */

#pragma once

#include <cmath>
#include <limits>

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"

#include "CubicHorizontalKernel.h"
#include "CubicVerticalKernel.h"

namespace atlas {
namespace interpolation {
namespace method {

class Cubic3DKernel {
public:
    Cubic3DKernel( const functionspace::StructuredColumns& fs, const util::Config& config = util::NoConfig() ) {
        src_ = fs;
        ATLAS_ASSERT( src_ );
        ATLAS_ASSERT( src_.halo() >= 2 );
        ATLAS_ASSERT( src_.vertical().size() );
        horizontal_interpolation_ = CubicHorizontalKernel( src_, config );
        vertical_interpolation_   = CubicVerticalKernel( fs.vertical(), config );
        limiter_                  = config.getBool( "limiter", false );
    }

private:
    functionspace::StructuredColumns src_;
    CubicHorizontalKernel horizontal_interpolation_;
    CubicVerticalKernel vertical_interpolation_;
    bool limiter_{false};

public:
    static std::string className() { return "Cubic3DKernel"; }
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

        if ( limiter_ ) {
            limit_scalar( output, index, stencil, input );
        }
        return output;
    }

    template <typename array_t, typename stencil_t>
    typename std::enable_if<( array_t::RANK == 2 ), void>::type limit_scalar(
        typename array_t::value_type& output, const std::array<std::array<idx_t, 4>, 4>& index,
        const stencil_t& stencil, const array_t& input ) const {
        using Scalar = typename array_t::value_type;
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x
        idx_t k = stencil.k_interval();
        idx_t k1, k2;
        if ( k < 0 ) {
            k1 = k2 = 0;
        }
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
        if ( output < minval ) {
            output = minval;
        }
        else if ( output > maxval ) {
            output = maxval;
        }
    }

    template <typename Value>
    struct OutputView1D {
        template <typename Int>
        Value& operator()( Int v ) {
            return data_[v];
        }
        template <typename Int>
        Value& operator[]( Int v ) {
            return data_[v];
        }
        static constexpr int RANK{1};
        OutputView1D( Value* data ) : data_( data ) {}
        using value_type = Value;

        Value* data_;
    };

    template <typename Value>
    OutputView1D<Value> make_outputview( Value* data ) const {
        return OutputView1D<Value>( data );
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<( InputArray::RANK == 3 ), void>::type interpolate_vars( const stencil_t& stencil,
                                                                                     const weights_t& weights,
                                                                                     const InputArray& input,
                                                                                     OutputArray& output,
                                                                                     const idx_t nvar ) const {
        using Value = typename InputArray::value_type;

        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& wj = weights.weights_j;
        const auto& wk = weights.weights_k;

        const Value* _input_;

        for ( idx_t v = 0; v < nvar; ++v ) {
            output[v] = 0.;
        }

        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& wi = weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                const idx_t n   = src_.index( stencil.i( i, j ), stencil.j( j ) );
                const Value wij = wi[i] * wj[j];
                for ( idx_t k = 0; k < stencil_width(); ++k ) {
                    const Value w  = wij * wk[k];
                    const idx_t kk = stencil.k( k );
                    _input_        = &( input( n, kk, 0 ) );  // Assumption that input.stride(2) == 1
                    for ( idx_t v = 0; v < nvar; ++v ) {
                        output[v] += w * _input_[v];
                    }
                }
                index[j][i] = n;
            }
        }

        if ( limiter_ ) {
            limit_vars( index, stencil, input, output, nvar );
        }
    }

    template <typename InputArray, typename OutputArray, typename stencil_t>
    typename std::enable_if<( InputArray::RANK == 3 ), void>::type limit_vars(
        const std::array<std::array<idx_t, 4>, 4>& index, const stencil_t& stencil, const InputArray& input,
        OutputArray& output, const idx_t nvar ) const {
        // Limit output to max/min of values in stencil marked by '*'
        //         x        x        x         x
        //              x     *-----*     x
        //                   /   P  |
        //          x       *------ *        x
        //        x        x        x         x

        using Value = typename InputArray::value_type;

        const idx_t k = stencil.k_interval();
        idx_t k1, k2;
        if ( k < 0 ) {
            k1 = k2 = stencil.k( 0 );
        }
        else if ( k > 2 ) {
            k1 = k2 = stencil.k( 3 );
        }
        else {
            k1 = stencil.k( k );
            k2 = k1 + 1;
        }

        for ( idx_t v = 0; v < nvar; ++v ) {
            Value limited = output[v];
            Value maxval  = std::numeric_limits<Value>::lowest();
            Value minval  = std::numeric_limits<Value>::max();
            for ( idx_t j = 1; j < 3; ++j ) {
                for ( idx_t i = 1; i < 3; ++i ) {
                    idx_t n = index[j][i];

                    Value f1 = input( n, k1, v );
                    Value f2 = input( n, k2, v );

                    maxval = std::max( maxval, f1 );
                    maxval = std::max( maxval, f2 );
                    minval = std::min( minval, f1 );
                    minval = std::min( minval, f2 );
                }
            }
            if ( limited < minval ) {
                limited = minval;
            }
            else if ( limited > maxval ) {
                limited = maxval;
            }
            output[v] = limited;
        }
    }


    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<( InputArray::RANK == 2 && OutputArray::RANK == 1 ), void>::type interpolate(
        const stencil_t& stencil, const weights_t& weights, const InputArray& input, OutputArray& output,
        idx_t r ) const {
        output( r ) = interpolate( stencil, weights, input );
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<( InputArray::RANK == 2 && OutputArray::RANK == 2 ), void>::type interpolate(
        const stencil_t& stencil, const weights_t& weights, const InputArray& input, OutputArray& output, idx_t r,
        idx_t k ) const {
        output( r, k ) = interpolate( stencil, weights, input );
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<( InputArray::RANK == 3 && OutputArray::RANK == 3 ), void>::type interpolate(
        const stencil_t& stencil, const weights_t& weights, const InputArray& input, OutputArray& output, idx_t r,
        idx_t k ) const {
        auto output_vars = make_outputview( &output( r, k, 0 ) );
        interpolate_vars( stencil, weights, input, output_vars, output.shape( 2 ) );
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<( InputArray::RANK == 2 && OutputArray::RANK == 3 ), void>::type interpolate(
        const stencil_t&, const weights_t&, const InputArray&, OutputArray&, idx_t /*r*/, idx_t /*k*/ ) const {
        ATLAS_NOTIMPLEMENTED;
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<( InputArray::RANK == 3 && OutputArray::RANK == 1 ), void>::type interpolate(
        const stencil_t&, const weights_t&, const InputArray&, OutputArray&, idx_t /*r*/ ) const {
        ATLAS_NOTIMPLEMENTED;
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<( InputArray::RANK == 3 && OutputArray::RANK == 1 ), void>::type interpolate(
        const stencil_t&, const weights_t&, const InputArray&, OutputArray&, idx_t /*r*/, idx_t /*k*/ ) const {
        ATLAS_NOTIMPLEMENTED;
    }

    template <typename stencil_t, typename weights_t, typename InputArray, typename OutputArray>
    typename std::enable_if<( InputArray::RANK == 3 && OutputArray::RANK == 2 ), void>::type interpolate(
        const stencil_t&, const weights_t&, const InputArray&, OutputArray&, idx_t /*r*/, idx_t /*k*/ ) const {
        ATLAS_NOTIMPLEMENTED;
    }
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
