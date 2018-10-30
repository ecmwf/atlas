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

#include "atlas/interpolation/method/Method.h"

#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/memory/NonCopyable.h"

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {

class CubicStructured2D : public Method {
public:
    CubicStructured2D( const Config& config ) : Method( config ), matrix_free_{false} {
        config.get( "matrix_free", matrix_free_ );
    }

    virtual ~CubicStructured2D() override {}

    virtual void setup( const FunctionSpace& source, const FunctionSpace& target ) override;

    virtual void print( std::ostream& ) const override;

protected:
    void setup( const FunctionSpace& source );

protected:
    Field target_xy_;
    Field target_ghost_;

    FunctionSpace source_;
    FunctionSpace target_;

    bool matrix_free_;

private:
    functionspace::StructuredColumns src_;
    ComputeHorizontalStencil compute_horizontal_stencil_;
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }
    bool limiter_{false};

public:
    using Stencil = HorizontalStencil<4>;
    struct Weights {
        std::array<std::array<double, 4>, 4> weights_i;
        std::array<double, 4> weights_j;
    };

public:
    template <typename stencil_t>
    void compute_stencil( const double x, const double y, stencil_t& stencil ) const {
        compute_horizontal_stencil_( x, y, stencil );
    }

    template <typename weights_t>
    void compute_weights( const double x, const double y, weights_t& weights ) const {
        HorizontalStencil<stencil_width()> stencil;
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
    void interpolate( const stencil_t& stencil, const weights_t& weights, const array_t& input, double& output ) {
        std::array<std::array<idx_t, stencil_width()>, stencil_width()> index;
        const auto& weights_j = weights.weights_j;
        output                = 0.;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& weights_i = weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t n = src_.index( stencil.i( i, j ), stencil.j( j ) );
                output += weights_i[i] * weights_j[j] * input[n];
                index[j][i] = n;
            }
        }

        if ( limiter_ ) { limit( output, index, input ); }
    }

    template <typename array_t>
    void limit( double& output, const std::array<std::array<idx_t, 4>, 4>& index, const array_t& input ) {
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
                idx_t n    = index[j][i];
                double val = input[n];
                maxval     = std::max( maxval, val );
                minval     = std::min( minval, val );
            }
        }
        output = std::min( maxval, std::max( minval, output ) );
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

    void insert_triplets( const idx_t row, const PointXY& p, Triplets& triplets, WorkSpace& ws ) {
        insert_triplets( row, p.x(), p.y(), triplets, ws );
    }

    void insert_triplets( const idx_t row, const double x, const double y, Triplets& triplets, WorkSpace& ws ) {
        compute_horizontal_stencil_( x, y, ws.stencil );
        compute_weights( x, y, ws.stencil, ws.weights );
        const auto& wj = ws.weights.weights_j;
        for ( idx_t j = 0; j < stencil_width(); ++j ) {
            const auto& wi = ws.weights.weights_i[j];
            for ( idx_t i = 0; i < stencil_width(); ++i ) {
                idx_t col = src_.index( ws.stencil.i( i, j ), ws.stencil.j( j ) );
                double w  = wi[i] * wj[j];
                triplets.emplace_back( row, col, w );
            }
        }
    }
};

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
