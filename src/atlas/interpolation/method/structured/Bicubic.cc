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

#include "atlas/interpolation/method/structured/Bicubic.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/geometry/Point3.h"
#include "eckit/log/Plural.h"
#include "eckit/log/ProgressTimer.h"
#include "eckit/log/Seconds.h"
#include "eckit/mpi/Comm.h"

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Stencil.h"
#include "atlas/grid/StencilComputer.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/Buffer.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace interpolation {
namespace method {


namespace {

MethodBuilder<Bicubic> __builder1( "structured-bicubic" );
MethodBuilder<Bicubic> __builder2( "bicubic" );

}  // namespace

namespace detail {

class BiCubicKernel {
    using Triplet  = eckit::linalg::Triplet;
    using Triplets = std::vector<Triplet>;

public:
    BiCubicKernel( const functionspace::StructuredColumns& fs ) {
        src_ = fs;
        ASSERT( src_ );
        ASSERT( src_.halo() >= 2 );
        compute_horizontal_stencil_ = ComputeHorizontalStencil( src_.grid(), stencil_width() );
    }

private:
    functionspace::StructuredColumns src_;
    ComputeHorizontalStencil compute_horizontal_stencil_;
    bool limiter_{false};
    static constexpr idx_t stencil_width() { return 4; }
    static constexpr idx_t stencil_size() { return stencil_width() * stencil_width(); }

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

}  // namespace detail

Bicubic::Bicubic( const Method::Config& config ) : Method( config ), matrix_free_{false} {
    config.get( "matrix_free", matrix_free_ );
}

void Bicubic::setup( const Grid& source, const Grid& target ) {
    if ( mpi::comm().size() > 1 ) { NOTIMP; }


    ASSERT( grid::StructuredGrid( source ) );
    FunctionSpace source_fs = functionspace::StructuredColumns( source, option::halo( 2 ) );
    FunctionSpace target_fs = functionspace::PointCloud( target );

    setup( source_fs, target_fs );
}

void Bicubic::setup( const FunctionSpace& source, const FunctionSpace& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::Bicubic::setup()" );

    source_ = source;
    target_ = target;

    if ( functionspace::NodeColumns tgt = target ) {
        target_lonlat_ = tgt.mesh().nodes().lonlat();
        target_ghost_  = tgt.mesh().nodes().ghost();
    }
    else if ( functionspace::PointCloud tgt = target ) {
        target_lonlat_ = tgt.lonlat();
        target_ghost_  = tgt.ghost();
    }
    else {
        NOTIMP;
    }

    setup( source );
}

void Bicubic::print( std::ostream& out ) const {
    ASSERT( not matrix_.empty() );

    functionspace::NodeColumns src( source_ );
    functionspace::NodeColumns tgt( target_ );
    if ( not tgt ) NOTIMP;
    auto gidx_src = array::make_view<gidx_t, 1>( src.nodes().global_index() );

    ASSERT( tgt.nodes().size() == idx_t( matrix_.rows() ) );


    auto field_stencil_points_loc  = tgt.createField<gidx_t>( option::variables( 16 ) );
    auto field_stencil_weigths_loc = tgt.createField<double>( option::variables( 16 ) );
    auto field_stencil_size_loc    = tgt.createField<int>();

    auto stencil_points_loc  = array::make_view<gidx_t, 2>( field_stencil_points_loc );
    auto stencil_weights_loc = array::make_view<double, 2>( field_stencil_weigths_loc );
    auto stencil_size_loc    = array::make_view<int, 1>( field_stencil_size_loc );
    stencil_size_loc.assign( 0 );

    for ( Matrix::const_iterator it = matrix_.begin(); it != matrix_.end(); ++it ) {
        idx_t p                     = idx_t( it.row() );
        idx_t& i                    = stencil_size_loc( p );
        stencil_points_loc( p, i )  = gidx_src( it.col() );
        stencil_weights_loc( p, i ) = *it;
        ++i;
    }


    gidx_t global_size = tgt.gather().glb_dof();

    auto field_stencil_points_glb  = tgt.createField<gidx_t>( option::variables( 16 ) | option::global( 0 ) );
    auto field_stencil_weights_glb = tgt.createField<double>( option::variables( 16 ) | option::global( 0 ) );
    auto field_stencil_size_glb    = tgt.createField<int>( option::global( 0 ) );


    auto stencil_points_glb  = array::make_view<gidx_t, 2>( field_stencil_points_glb );
    auto stencil_weights_glb = array::make_view<double, 2>( field_stencil_weights_glb );
    auto stencil_size_glb    = array::make_view<int, 1>( field_stencil_size_glb );

    tgt.gather().gather( stencil_size_loc, stencil_size_glb );
    tgt.gather().gather( stencil_points_loc, stencil_points_glb );
    tgt.gather().gather( stencil_weights_loc, stencil_weights_glb );

    if ( mpi::comm().rank() == 0 ) {
        int precision = std::numeric_limits<double>::max_digits10;
        for ( idx_t i = 0; i < global_size; ++i ) {
            out << std::setw( 10 ) << i + 1 << " : ";
            for ( idx_t j = 0; j < stencil_size_glb( i ); ++j ) {
                out << std::setw( 10 ) << stencil_points_glb( i, j );
            }
            for ( idx_t j = stencil_size_glb( i ); j < 16; ++j ) {
                out << "          ";
            }
            for ( idx_t j = 0; j < stencil_size_glb( i ); ++j ) {
                out << std::setw( precision + 5 ) << std::left << std::setprecision( precision )
                    << stencil_weights_glb( i, j );
            }
            out << std::endl;
        }
    }
}

void Bicubic::setup( const FunctionSpace& source ) {
    kernel_.reset( new Kernel( source ) );

    if ( not matrix_free_ ) {
        idx_t inp_npts = source.size();
        idx_t out_npts = target_lonlat_.shape( 0 );

        auto ghost  = array::make_view<int, 1>( target_ghost_ );
        auto lonlat = array::make_view<double, 2>( target_lonlat_ );

        auto triplets = kernel_->allocate_triplets( out_npts );

        constexpr NormaliseLongitude normalise;
        ATLAS_TRACE_SCOPE( "Precomputing interpolation matrix" ) {
            atlas_omp_parallel {
                Kernel::WorkSpace workspace;
                atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
                    if ( not ghost( n ) ) {
                        PointLonLat p{lonlat( n, LON ), lonlat( n, LAT )};
                        normalise( p );
                        kernel_->insert_triplets( n, p, triplets, workspace );
                    }
                }
            }
            // fill sparse matrix and return
            Matrix A( out_npts, inp_npts, triplets );
            matrix_.swap( A );
        }
    }
}

void Bicubic::execute( const Field& src_field, Field& tgt_field ) const {
    if ( not matrix_free_ ) {
        Method::execute( src_field, tgt_field );
        return;
    }

    if ( src_field.dirty() ) { source().haloExchange( const_cast<Field&>( src_field ) ); }

    ATLAS_TRACE( "atlas::interpolation::method::Bicubic::execute()" );

    array::DataType datatype = src_field.datatype();
    int rank                 = src_field.rank();

    ASSERT( tgt_field.datatype() == datatype );
    ASSERT( tgt_field.rank() == rank );

    if ( datatype.kind() == array::DataType::KIND_REAL64 && rank == 1 ) {
        execute_impl<double, 1>( src_field, tgt_field );
    }
    if ( datatype.kind() == array::DataType::KIND_REAL32 && rank == 1 ) {
        execute_impl<float, 1>( src_field, tgt_field );
    }
    if ( datatype.kind() == array::DataType::KIND_REAL64 && rank == 2 ) {
        execute_impl<double, 2>( src_field, tgt_field );
    }
    if ( datatype.kind() == array::DataType::KIND_REAL32 && rank == 2 ) {
        execute_impl<float, 2>( src_field, tgt_field );
    }

    tgt_field.set_dirty();
}

void Bicubic::execute( const FieldSet& src_fields, FieldSet& tgt_fields ) const {
    if ( not matrix_free_ ) {
        Method::execute( src_fields, tgt_fields );
        return;
    }

    ATLAS_TRACE( "atlas::interpolation::method::Bicubic::execute()" );

    const idx_t N = src_fields.size();
    ASSERT( N == tgt_fields.size() );

    if ( N == 0 ) return;

    for ( idx_t i = 0; i < N; ++i ) {
        if ( src_fields[i].dirty() ) { source().haloExchange( const_cast<Field&>( src_fields[i] ) ); }
    }

    array::DataType datatype = src_fields[0].datatype();
    int rank                 = src_fields[0].rank();

    for ( idx_t i = 0; i < N; ++i ) {
        ASSERT( src_fields[i].datatype() == datatype );
        ASSERT( src_fields[i].rank() == rank );
        ASSERT( tgt_fields[i].datatype() == datatype );
        ASSERT( tgt_fields[i].rank() == rank );
    }

    if ( datatype.kind() == array::DataType::KIND_REAL64 && rank == 1 ) {
        execute_impl<double, 1>( src_fields, tgt_fields );
    }
    if ( datatype.kind() == array::DataType::KIND_REAL32 && rank == 1 ) {
        execute_impl<float, 1>( src_fields, tgt_fields );
    }
    if ( datatype.kind() == array::DataType::KIND_REAL64 && rank == 2 ) {
        execute_impl<double, 2>( src_fields, tgt_fields );
    }
    if ( datatype.kind() == array::DataType::KIND_REAL32 && rank == 2 ) {
        execute_impl<float, 2>( src_fields, tgt_fields );
    }

    for ( idx_t i = 0; i < N; ++i ) {
        tgt_fields[i].set_dirty();
    }
}

template <typename Value, int Rank>
void Bicubic::execute_impl( const FieldSet& src_fields, FieldSet& tgt_fields ) const {
    const idx_t N  = src_fields.size();
    idx_t out_npts = target_lonlat_.shape( 0 );

    auto ghost  = array::make_view<int, 1>( target_ghost_ );
    auto lonlat = array::make_view<double, 2>( target_lonlat_ );

    std::vector<array::ArrayView<Value, Rank> > src_view;
    std::vector<array::ArrayView<Value, Rank> > tgt_view;
    src_view.reserve( N );
    tgt_view.reserve( N );

    for ( idx_t i = 0; i < N; ++i ) {
        src_view.emplace_back( array::make_view<Value, Rank>( src_fields[i] ) );
        tgt_view.emplace_back( array::make_view<Value, Rank>( tgt_fields[i] ) );
    }

    constexpr NormaliseLongitude normalise( 0., 360. );  // includes 360 as well!
    atlas_omp_parallel {
        Kernel::Stencil stencil;
        Kernel::Weights weights;
        atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
            if ( not ghost( n ) ) {
                PointLonLat p{lonlat( n, LON ), lonlat( n, LAT )};
                normalise( p );
                kernel_->compute_stencil( p.lon(), p.lat(), stencil );
                kernel_->compute_weights( p.lon(), p.lat(), stencil, weights );
                for ( idx_t i = 0; i < N; ++i ) {
                    kernel_->interpolate( stencil, weights, src_view[i], tgt_view[i], n );
                }
            }
        }
    }
}

template <typename Value, int Rank>
void Bicubic::execute_impl( const Field& src_field, Field& tgt_field ) const {
    idx_t out_npts = target_lonlat_.shape( 0 );

    auto ghost  = array::make_view<int, 1>( target_ghost_ );
    auto lonlat = array::make_view<double, 2>( target_lonlat_ );

    const auto src_view = array::make_view<Value, Rank>( src_field );
    auto tgt_view       = array::make_view<Value, Rank>( tgt_field );

    constexpr NormaliseLongitude normalise( 0., 360. );  // includes 360 as well!
    atlas_omp_parallel {
        Kernel::Stencil stencil;
        Kernel::Weights weights;
        atlas_omp_for( idx_t n = 0; n < out_npts; ++n ) {
            if ( not ghost( n ) ) {
                PointLonLat p{lonlat( n, LON ), lonlat( n, LAT )};
                normalise( p );
                kernel_->compute_stencil( p.lon(), p.lat(), stencil );
                kernel_->compute_weights( p.lon(), p.lat(), stencil, weights );
                kernel_->interpolate( stencil, weights, src_view, tgt_view, n );
            }
        }
    }
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
