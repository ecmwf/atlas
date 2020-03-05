/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/interpolation/method/Method.h"

#include "eckit/linalg/LinearAlgebra.h"
#include "eckit/linalg/Vector.h"
#include "eckit/log/Timer.h"
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/thread/Once.h"

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

namespace atlas {
namespace interpolation {

void Method::check_compatibility( const Field& src, const Field& tgt ) const {
    ATLAS_ASSERT( src.datatype() == tgt.datatype() );
    ATLAS_ASSERT( src.rank() == tgt.rank() );
    ATLAS_ASSERT( src.levels() == tgt.levels() );
    ATLAS_ASSERT( src.variables() == tgt.variables() );

    ATLAS_ASSERT( !matrix_.empty() );
    ATLAS_ASSERT( tgt.shape( 0 ) >= static_cast<idx_t>( matrix_.rows() ) );
    ATLAS_ASSERT( src.shape( 0 ) >= static_cast<idx_t>( matrix_.cols() ) );
}

template <typename Value>
void Method::interpolate_field( const Field& src, Field& tgt ) const {
    check_compatibility( src, tgt );
    if ( src.rank() == 1 ) {
        interpolate_field_rank1<Value>( src, tgt );
    }
    if ( src.rank() == 2 ) {
        interpolate_field_rank2<Value>( src, tgt );
    }
    if ( src.rank() == 3 ) {
        interpolate_field_rank3<Value>( src, tgt );
    }
}

template <typename Value>
void Method::interpolate_field_rank1( const Field& src, Field& tgt ) const {
    const auto outer  = matrix_.outer();
    const auto index  = matrix_.inner();
    const auto weight = matrix_.data();
    idx_t rows        = static_cast<idx_t>( matrix_.rows() );

    if ( use_eckit_linalg_spmv_ ) {
        if ( src.datatype() != array::make_datatype<double>() ) {
            throw_NotImplemented( "Only double precision interpolation is currently implemented with eckit backend",
                                  Here() );
        }
        ATLAS_ASSERT( src.contiguous() );
        ATLAS_ASSERT( tgt.contiguous() );

        eckit::linalg::Vector v_src( const_cast<double*>( array::make_view<double, 1>( src ).data() ), src.shape( 0 ) );
        eckit::linalg::Vector v_tgt( array::make_view<double, 1>( tgt ).data(), tgt.shape( 0 ) );
        eckit::linalg::LinearAlgebra::backend().spmv( matrix_, v_src, v_tgt );
    }
    else {
        auto v_src = array::make_view<Value, 1>( src );
        auto v_tgt = array::make_view<Value, 1>( tgt );

        atlas_omp_parallel_for( idx_t r = 0; r < rows; ++r ) {
            v_tgt( r ) = 0.;
            for ( idx_t c = outer[r]; c < outer[r + 1]; ++c ) {
                idx_t n = index[c];
                Value w = static_cast<Value>( weight[c] );
                v_tgt( r ) += w * v_src( n );
            }
        }
    }
}

template <typename Value>
void Method::interpolate_field_rank2( const Field& src, Field& tgt ) const {
    const auto outer  = matrix_.outer();
    const auto index  = matrix_.inner();
    const auto weight = matrix_.data();
    idx_t rows        = static_cast<idx_t>( matrix_.rows() );

    auto v_src = array::make_view<Value, 2>( src );
    auto v_tgt = array::make_view<Value, 2>( tgt );

    idx_t Nk = src.shape( 1 );

    atlas_omp_parallel_for( idx_t r = 0; r < rows; ++r ) {
        for ( idx_t k = 0; k < Nk; ++k ) {
            v_tgt( r, k ) = 0.;
        }
        for ( idx_t c = outer[r]; c < outer[r + 1]; ++c ) {
            idx_t n = index[c];
            Value w = static_cast<Value>( weight[c] );
            for ( idx_t k = 0; k < Nk; ++k ) {
                v_tgt( r, k ) += w * v_src( n, k );
            }
        }
    }
}


template <typename Value>
void Method::interpolate_field_rank3( const Field& src, Field& tgt ) const {
    const auto outer  = matrix_.outer();
    const auto index  = matrix_.inner();
    const auto weight = matrix_.data();
    idx_t rows        = static_cast<idx_t>( matrix_.rows() );

    auto v_src = array::make_view<Value, 3>( src );
    auto v_tgt = array::make_view<Value, 3>( tgt );

    idx_t Nk = src.shape( 1 );
    idx_t Nl = src.shape( 2 );

    atlas_omp_parallel_for( idx_t r = 0; r < rows; ++r ) {
        for ( idx_t k = 0; k < Nk; ++k ) {
            for ( idx_t l = 0; l < Nl; ++l ) {
                v_tgt( r, k, l ) = 0.;
            }
        }
        for ( idx_t c = outer[r]; c < outer[r + 1]; ++c ) {
            idx_t n = index[c];
            Value w = static_cast<Value>( weight[c] );
            for ( idx_t k = 0; k < Nk; ++k ) {
                for ( idx_t l = 0; l < Nl; ++l ) {
                    v_tgt( r, k, l ) += w * v_src( n, k, l );
                }
            }
        }
    }
}

Method::Method( const Method::Config& config ) {
    std::string spmv = "";
    config.get( "spmv", spmv );
    use_eckit_linalg_spmv_ = ( spmv == "eckit" );
}

void Method::setup( const FunctionSpace& source, const FunctionSpace& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::Method::setup(FunctionSpace, FunctionSpace)" );
    this->do_setup( source, target );
}

void Method::setup( const Grid& source, const Grid& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::Method::setup(Grid, Grid)" );
    this->do_setup( source, target );
}

void Method::setup( const FunctionSpace& source, const Field& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::Method::setup(FunctionSpace, Field)" );
    this->do_setup( source, target );
}

void Method::setup( const FunctionSpace& source, const FieldSet& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::Method::setup(FunctionSpace, FieldSet)" );
    this->do_setup( source, target );
}

void Method::execute( const FieldSet& source, FieldSet& target ) const {
    ATLAS_TRACE( "atlas::interpolation::method::Method::execute(FieldSet, FieldSet)" );
    this->do_execute( source, target );
}

void Method::execute( const Field& source, Field& target ) const {
    ATLAS_TRACE( "atlas::interpolation::method::Method::execute(Field, Field)" );
    this->do_execute( source, target );
}

void Method::do_setup( const FunctionSpace& /*source*/, const Field& /*target*/ ) {
    ATLAS_NOTIMPLEMENTED;
}

void Method::do_setup( const FunctionSpace& /*source*/, const FieldSet& /*target*/ ) {
    ATLAS_NOTIMPLEMENTED;
}

void Method::do_execute( const FieldSet& fieldsSource, FieldSet& fieldsTarget ) const {
    ATLAS_TRACE( "atlas::interpolation::method::Method::do_execute()" );

    const idx_t N = fieldsSource.size();
    ATLAS_ASSERT( N == fieldsTarget.size() );

    for ( idx_t i = 0; i < fieldsSource.size(); ++i ) {
        Log::debug() << "Method::do_execute() on field " << ( i + 1 ) << '/' << N << "..." << std::endl;
        Method::do_execute( fieldsSource[i], fieldsTarget[i] );
    }
}

void Method::do_execute( const Field& src, Field& tgt ) const {
    haloExchange( src );

    ATLAS_TRACE( "atlas::interpolation::method::Method::do_execute()" );

    if ( src.datatype().kind() == array::DataType::KIND_REAL64 ) {
        interpolate_field<double>( src, tgt );
    }
    if ( src.datatype().kind() == array::DataType::KIND_REAL32 ) {
        interpolate_field<float>( src, tgt );
    }

    tgt.set_dirty();
}

void Method::normalise( Triplets& triplets ) {
    // sum all calculated weights for normalisation
    double sum = 0.0;

    for ( size_t j = 0; j < triplets.size(); ++j ) {
        sum += triplets[j].value();
    }

    // now normalise all weights according to the total
    const double invSum = 1.0 / sum;
    for ( size_t j = 0; j < triplets.size(); ++j ) {
        triplets[j].value() *= invSum;
    }
}

void Method::haloExchange( const FieldSet& fields ) const {
    for ( auto& field : fields ) {
        haloExchange( field );
    }
}
void Method::haloExchange( const Field& field ) const {
    if ( field.dirty() ) {
        source().haloExchange( field );
    }
}

}  // namespace interpolation
}  // namespace atlas
