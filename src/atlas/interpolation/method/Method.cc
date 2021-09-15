/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <memory>

#include "atlas/interpolation/method/Method.h"


#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/linalg/sparse.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

using namespace atlas::linalg;
namespace atlas {
namespace interpolation {

namespace {

template <typename Value>
void set_missing_values_rank1( Field& tgt, const std::vector<idx_t>& missing, const Value& missing_value ) {
    auto tgt_v = array::make_view<Value, 1>( tgt );
    for ( auto i : missing ) {
        tgt_v( i ) = missing_value;
    }
}

template <typename Value>
void set_missing_values_rank2( Field& tgt, const std::vector<idx_t>& missing, const Value& missing_value ) {
    auto tgt_v     = array::make_view<Value, 2>( tgt );
    const idx_t Nj = tgt_v.shape( 1 );
    for ( auto i : missing ) {
        for ( idx_t j = 0; j < Nj; ++j ) {
            tgt_v( i, j ) = missing_value;
        }
    }
}

template <typename Value>
void set_missing_values_rank3( Field& tgt, const std::vector<idx_t>& missing, const Value& missing_value ) {
    auto tgt_v     = array::make_view<Value, 3>( tgt );
    const idx_t Nj = tgt_v.shape( 1 );
    const idx_t Nk = tgt_v.shape( 2 );
    for ( auto i : missing ) {
        for ( idx_t j = 0; j < Nj; ++j ) {
            for ( idx_t k = 0; k < Nk; ++k ) {
                tgt_v( i, j, k ) = missing_value;
            }
        }
    }
}

template <typename Value>
void set_missing_values_T( Field& tgt, const std::vector<idx_t>& missing ) {
    Value missing_value = tgt.metadata().get<Value>( "missing_value" );
    if ( tgt.rank() == 1 ) {
        set_missing_values_rank1( tgt, missing, missing_value );
    }
    else if ( tgt.rank() == 2 ) {
        set_missing_values_rank2( tgt, missing, missing_value );
    }
    else if ( tgt.rank() == 3 ) {
        set_missing_values_rank3( tgt, missing, missing_value );
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void set_missing_values( Field& tgt, const std::vector<idx_t>& missing ) {
    if ( missing.empty() ) {
        return;
    }
    if ( tgt.datatype().kind() == array::DataType::KIND_REAL64 ) {
        set_missing_values_T<double>( tgt, missing );
    }
    else if ( tgt.datatype().kind() == array::DataType::KIND_REAL32 ) {
        set_missing_values_T<float>( tgt, missing );
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

}  // anonymous namespace

template <typename Value>
void Method::interpolate_field_rank1( const Field& src, Field& tgt, const Matrix& W ) const {
    if ( use_eckit_linalg_spmv_ ) {
        if ( src.datatype() != array::make_datatype<double>() ) {
            throw_NotImplemented( "Only double precision interpolation is currently implemented with eckit backend",
                                  Here() );
        }
        auto src_v = array::make_view<double, 1>( src );
        auto tgt_v = array::make_view<double, 1>( tgt );
        sparse_matrix_multiply( W, src_v, tgt_v, sparse::backend::eckit_linalg() );
    }
    else {
        auto src_v = array::make_view<Value, 1>( src );
        auto tgt_v = array::make_view<Value, 1>( tgt );
        sparse_matrix_multiply( W, src_v, tgt_v, sparse::backend::omp() );
    }
}

template <typename Value>
void Method::interpolate_field_rank2( const Field& src, Field& tgt, const Matrix& W ) const {
    auto src_v = array::make_view<Value, 2>( src );
    auto tgt_v = array::make_view<Value, 2>( tgt );
    sparse_matrix_multiply( W, src_v, tgt_v, sparse::backend::omp() );
}


template <typename Value>
void Method::interpolate_field_rank3( const Field& src, Field& tgt, const Matrix& W ) const {
    auto src_v = array::make_view<Value, 3>( src );
    auto tgt_v = array::make_view<Value, 3>( tgt );
    sparse_matrix_multiply( W, src_v, tgt_v, sparse::backend::omp() );
}

template <typename Value>
void Method::adjoint_interpolate_field_rank1( Field& src, Field& tgt, const Matrix& W ) const {

    Field tmp( src.name(), src.datatype().kind(), src.shape() );
    tmp.set_levels( src.levels() );

    auto tmp_v = array::make_view<double, 1>( tmp );
    auto src_v = array::make_view<double, 1>( src );
    auto tgt_v = array::make_view<double, 1>( tgt );

    for ( idx_t t = 0; t < tmp.shape(0); ++t ) {
        tmp_v(t) = 0.;
    }

    if ( use_eckit_linalg_spmv_ ) {
        if ( src.datatype() != array::make_datatype<double>() ) {
            throw_NotImplemented( "Only double precision interpolation is currently implemented with eckit backend",
                                  Here() );
        }
        sparse_matrix_multiply( W, tgt_v, tmp_v, sparse::backend::eckit_linalg() );
    }
    else {
        sparse_matrix_multiply( W, tgt_v, tmp_v, sparse::backend::omp() );
    }

    for ( idx_t t = 0; t < tmp.shape(0); ++t ) {
        src_v(t) += tmp_v(t);
        tgt_v(t) = 0.;
    }
}

template <typename Value>
void Method::adjoint_interpolate_field_rank2( Field& src, Field& tgt, const Matrix& W ) const {
    Field tmp( src.name(), src.datatype().kind(), src.shape() );
    tmp.set_levels( src.levels() );

    auto tmp_v = array::make_view<Value, 2>( tmp );
    auto src_v = array::make_view<Value, 2>( src );
    auto tgt_v = array::make_view<Value, 2>( tgt );

    for ( idx_t t = 0; t < tmp.shape(0); ++t ) {
        for ( idx_t k = 0; k < tmp.shape(1); ++k ) {
            tmp_v(t, k) = 0.;
        }
    }

    sparse_matrix_multiply( W, tgt_v, tmp_v, sparse::backend::omp() );

    for ( idx_t t = 0; t < tmp.shape(0); ++t ) {
        for ( idx_t k = 0; k < tmp.shape(1); ++k ) {
            src_v(t, k) += tmp_v(t, k);
            tgt_v(t, k) = 0.;
        }
    }
}

template <typename Value>
void Method::adjoint_interpolate_field_rank3( Field& src, Field& tgt, const Matrix& W ) const {

    Field tmp( src.name(), src.datatype().kind(), src.shape() );
    tmp.set_levels( src.levels() );

    auto tmp_v = array::make_view<Value, 3>( tmp );
    auto src_v = array::make_view<Value, 3>( src );
    auto tgt_v = array::make_view<Value, 3>( tgt );

    for ( idx_t t = 0; t < tmp.shape(0); ++t ) {
        for ( idx_t j = 0; j < tmp.shape(1); ++j ) {
            for ( idx_t k = 0; k < tmp.shape(2); ++k ) {
                tmp_v(t, j, k) = 0.;
            }
        }
    }

    sparse_matrix_multiply( W, tgt_v, tmp_v, sparse::backend::omp() );

    for ( idx_t t = 0; t < tmp.shape(0); ++t ) {
        for ( idx_t j = 0; j < tmp.shape(1); ++j ) {
            for ( idx_t k = 0; k < tmp.shape(2); ++k ) {
                src_v(t, j, k) += tmp_v(t, j, k);
                tgt_v(t, j, k) = 0.;
            }
        }
    }
}

void Method::check_compatibility( const Field& src, const Field& tgt, const Matrix& W ) const {
    ATLAS_ASSERT( src.datatype() == tgt.datatype() );
    ATLAS_ASSERT( src.rank() == tgt.rank() );
    ATLAS_ASSERT( src.levels() == tgt.levels() );
    ATLAS_ASSERT( src.variables() == tgt.variables() );

    ATLAS_ASSERT( !W.empty() );
    ATLAS_ASSERT( tgt.shape( 0 ) >= static_cast<idx_t>( W.rows() ) );
    ATLAS_ASSERT( src.shape( 0 ) >= static_cast<idx_t>( W.cols() ) );
}

template <typename Value>
void Method::interpolate_field( const Field& src, Field& tgt, const Matrix& W ) const {
    check_compatibility( src, tgt, W );

    if ( src.rank() == 1 ) {
        interpolate_field_rank1<Value>( src, tgt, W );
    }
    else if ( src.rank() == 2 ) {
        interpolate_field_rank2<Value>( src, tgt, W );
    }
    else if ( src.rank() == 3 ) {
        interpolate_field_rank3<Value>( src, tgt, W );
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

template <typename Value>
void Method::adjoint_interpolate_field( Field& src, Field& tgt, const Matrix& W ) const {
    check_compatibility( tgt, src, W );

    if ( src.rank() == 1 ) {
        adjoint_interpolate_field_rank1<Value>( src, tgt, W );
    }
    else if ( src.rank() == 2 ) {
        adjoint_interpolate_field_rank2<Value>( src, tgt, W );
    }
    else if ( src.rank() == 3 ) {
        adjoint_interpolate_field_rank3<Value>( src, tgt, W );
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

Method::Method( const Method::Config& config ) {
    std::string spmv = "";
    config.get( "spmv", spmv );
    use_eckit_linalg_spmv_ = ( spmv == "eckit" );

    std::string non_linear;
    if ( config.get( "non_linear", non_linear ) ) {
        nonLinear_ = NonLinear( non_linear, config );
    }
    matrix_shared_ = std::make_shared<Matrix>();
    matrix_cache_  = interpolation::MatrixCache( matrix_shared_ );
    matrix_        = matrix_shared_.get();

    config.get("adjoint coefficients", create_adjoint_coeffs_);
}

void Method::setup( const FunctionSpace& source, const FunctionSpace& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::Method::setup(FunctionSpace, FunctionSpace)" );
    this->do_setup( source, target );

    if (create_adjoint_coeffs_) {
        Matrix tmp(*matrix_);

        // if interpolation is matrix free then matrix->nonZeros() will be zero.
        if (tmp.nonZeros() > 0) {
            matrix_transpose_ = tmp.transpose();
        }
    }
}

void Method::setup( const Grid& source, const Grid& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::Method::setup(Grid, Grid)" );
    this->do_setup( source, target, Cache() );
}

void Method::setup( const FunctionSpace& source, const Field& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::Method::setup(FunctionSpace, Field)" );
    this->do_setup( source, target );
}

void Method::setup( const FunctionSpace& source, const FieldSet& target ) {
    ATLAS_TRACE( "atlas::interpolation::method::Method::setup(FunctionSpace, FieldSet)" );
    this->do_setup( source, target );
}

void Method::setup( const Grid& source, const Grid& target, const Cache& cache ) {
    ATLAS_TRACE( "atlas::interpolation::method::Method::setup(Grid, Grid, Cache)" );
    this->do_setup( source, target, cache );
}

void Method::execute( const FieldSet& source, FieldSet& target ) const {
    ATLAS_TRACE( "atlas::interpolation::method::Method::execute(FieldSet, FieldSet)" );
    this->do_execute( source, target );
}

void Method::execute( const Field& source, Field& target ) const {
    ATLAS_TRACE( "atlas::interpolation::method::Method::execute(Field, Field)" );
    this->do_execute( source, target );
}

void Method::execute_adjoint( FieldSet& source, FieldSet& target ) const {
    ATLAS_TRACE( "atlas::interpolation::method::Method::execute(FieldSet, FieldSet)" );
    this->do_execute_adjoint( source, target );
}

void Method::execute_adjoint( Field& source, Field& target ) const {
    ATLAS_TRACE( "atlas::interpolation::method::Method::execute(Field, Field)" );
    this->do_execute_adjoint( source, target );
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
        Method::do_execute( fieldsSource[i], fieldsTarget[i] );
    }
}

void Method::do_execute( const Field& src, Field& tgt ) const {
    ATLAS_TRACE( "atlas::interpolation::method::Method::do_execute()" );

    haloExchange( src );

    // non-linearities: a non-empty M matrix contains the corrections applied to matrix_
    Matrix M;
    if ( not matrix_->empty() && nonLinear_( src ) ) {
        Matrix W( *matrix_ );  // copy (a big penalty -- copy-on-write would definitely be better)
        if ( nonLinear_->execute( W, src ) ) {
            M.swap( W );
        }
    }

    if ( src.datatype().kind() == array::DataType::KIND_REAL64 ) {
        interpolate_field<double>( src, tgt, M.empty() ? *matrix_ : M );
    }
    else if ( src.datatype().kind() == array::DataType::KIND_REAL32 ) {
        interpolate_field<float>( src, tgt, M.empty() ? *matrix_ : M );
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }

    // carry over missing value metadata
    if ( not tgt.metadata().has( "missing_value" ) ) {
        field::MissingValue mv_src( src );
        if ( mv_src ) {
            mv_src.metadata( tgt );
            ATLAS_ASSERT( field::MissingValue( tgt ) );
        }
        else if ( not missing_.empty() ) {
            if ( not tgt.metadata().has( "missing_value" ) ) {
                tgt.metadata().set( "missing_value", 9999. );
            }
            tgt.metadata().set( "missing_value_type", "equals" );
        }
    }

    // set missing values
    set_missing_values( tgt, missing_ );

    tgt.set_dirty();
}

void Method::do_execute_adjoint( FieldSet& fieldsSource, FieldSet& fieldsTarget ) const {
    ATLAS_TRACE( "atlas::interpolation::method::Method::do_execute_adjoint()" );

    const idx_t N = fieldsSource.size();
    ATLAS_ASSERT( N == fieldsTarget.size() );

    for ( idx_t i = 0; i < fieldsSource.size(); ++i ) {
        Method::do_execute_adjoint( fieldsSource[i], fieldsTarget[i] );
    }
}

void Method::do_execute_adjoint(Field& src, Field& tgt ) const {
    ATLAS_TRACE( "atlas::interpolation::method::Method::do_execute_adjoint()" );

    if ( nonLinear_( src ) ) {
        throw_NotImplemented( "Adjoint interpolation only works for interpolation schemes that are linear",
                              Here() );
    }

    if ( not missing_.empty() ) {
        throw_NotImplemented( "Adjoint Interpolation does not work for fields that have missing data. ",
                              Here() );
    }

    if ( matrix_transpose_.empty() ) {
        throw_AssertionFailed("Need to set 'adjoint coefficients' to true in config for adjoint interpolation to work");
    }

    if ( src.datatype().kind() == array::DataType::KIND_REAL64 ) {

        adjoint_interpolate_field<double>( src, tgt, matrix_transpose_ );
    }
    else if ( src.datatype().kind() == array::DataType::KIND_REAL32 ) {

        adjoint_interpolate_field<float>( src, tgt, matrix_transpose_ );
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }

    src.set_dirty();

    adjointHaloExchange( src );

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
    if ( field.dirty() && allow_halo_exchange_ ) {
        source().haloExchange( field );
    }
}

void Method::adjointHaloExchange( const FieldSet& fields ) const {
    for ( auto& field : fields ) {
        adjointHaloExchange( field );
    }
}
void Method::adjointHaloExchange( const Field& field ) const {
    if ( field.dirty() && allow_halo_exchange_ ) {
        source().adjointHaloExchange( field );
    }
}

interpolation::Cache Method::createCache() const {
    return matrix_cache_;
}


}  // namespace interpolation
}  // namespace atlas
