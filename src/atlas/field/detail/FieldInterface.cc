/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "atlas/field/Field.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/runtime/ErrorHandling.h"

namespace atlas {
namespace field {

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

namespace {
template <typename Value>
void atlas__Field__host_data_specf( FieldImpl* This, Value*& data, int& rank, int*& shapef, int*& stridesf ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); if ( This->datatype() != array::make_datatype<Value>() ) {
        throw eckit::Exception( "Datatype mismatch for accessing field data" );
    } This->array()
                                              .accMap();
                          data = This->host_data<Value>(); shapef = const_cast<int*>( This->shapef().data() );
                          stridesf = const_cast<int*>( This->stridesf().data() ); rank = This->shapef().size(); );
}
}  // namespace

extern "C" {

FieldImpl* atlas__Field__wrap_int_specf( const char* name, int data[], int rank, int shapef[], int stridesf[] ) {
    ATLAS_ERROR_HANDLING( array::ArrayShape shape; shape.resize( rank ); array::ArrayStrides strides;
                          strides.resize( rank ); idx_t jf = rank - 1; for ( int j = 0; j < rank; ++j ) {
                              shape[j]   = shapef[jf];
                              strides[j] = stridesf[jf];
                              --jf;
                          } FieldImpl * field;
                          {
                              Field wrapped( std::string( name ), data, array::ArraySpec( shape, strides ) );
                              field = wrapped.get();
                              field->attach();
                          } field->detach();
                          ASSERT( field ); return field; );
    return nullptr;
}

FieldImpl* atlas__Field__wrap_long_specf( const char* name, long data[], int rank, int shapef[], int stridesf[] ) {
    ATLAS_ERROR_HANDLING( array::ArrayShape shape; shape.resize( rank ); array::ArrayStrides strides;
                          strides.resize( rank ); idx_t jf = rank - 1; for ( int j = 0; j < rank; ++j ) {
                              shape[j]   = shapef[jf];
                              strides[j] = stridesf[jf];
                              --jf;
                          } FieldImpl * field;
                          {
                              Field wrapped( std::string( name ), data, array::ArraySpec( shape, strides ) );
                              field = wrapped.get();
                              field->attach();
                          } field->detach();
                          ASSERT( field ); return field; );
    return nullptr;
}

FieldImpl* atlas__Field__wrap_float_specf( const char* name, float data[], int rank, int shapef[], int stridesf[] ) {
    ATLAS_ERROR_HANDLING( array::ArrayShape shape; shape.resize( rank ); array::ArrayStrides strides;
                          strides.resize( rank ); idx_t jf = rank - 1; for ( int j = 0; j < rank; ++j ) {
                              shape[j]   = shapef[jf];
                              strides[j] = stridesf[jf];
                              --jf;
                          } FieldImpl * field;
                          {
                              Field wrapped( std::string( name ), data, array::ArraySpec( shape, strides ) );
                              field = wrapped.get();
                              field->attach();
                          } field->detach();
                          ASSERT( field ); return field; );
    return 0;
}

FieldImpl* atlas__Field__wrap_double_specf( const char* name, double data[], int rank, int shapef[], int stridesf[] ) {
    ATLAS_ERROR_HANDLING( array::ArrayShape shape; shape.resize( rank ); array::ArrayStrides strides;
                          strides.resize( rank ); idx_t jf = rank - 1; for ( int j = 0; j < rank; ++j ) {
                              shape[j]   = shapef[jf];
                              strides[j] = stridesf[jf];
                              --jf;
                          } FieldImpl * field;
                          {
                              Field wrapped( std::string( name ), data, array::ArraySpec( shape, strides ) );
                              field = wrapped.get();
                              field->attach();
                          } field->detach();
                          ASSERT( field ); return field; );
    return nullptr;
}

FieldImpl* atlas__Field__create( eckit::Parametrisation* params ) {
    ATLAS_ERROR_HANDLING( ASSERT( params ); FieldImpl * field; {
        Field f( *params );
        field = f.get();
        field->attach();
    } field->detach();

                          ASSERT( field ); return field; );
    return nullptr;
}

void atlas__Field__delete( FieldImpl* This ) {
    delete This;
}

const char* atlas__Field__name( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->name().c_str(); );
    return nullptr;
}

void atlas__Field__datatype( FieldImpl* This, char*& datatype, int& size, int& allocated ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); std::string s = This->datatype().str(); size = s.size() + 1;
                          datatype = new char[size]; strcpy( datatype, s.c_str() ); allocated = true; );
}

int atlas__Field__size( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->size(); );
    return 0;
}

int atlas__Field__rank( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->rank(); );
    return 0;
}

int atlas__Field__kind( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->datatype().kind(); );
    return 0;
}

double atlas__Field__bytes( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->bytes(); );
    return 0;
}

int atlas__Field__levels( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->levels(); );
    return 0;
}

util::Metadata* atlas__Field__metadata( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return &This->metadata(); );
    return nullptr;
}

int atlas__Field__has_functionspace( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return ( This->functionspace() != 0 ); );
    return 0;
}

const functionspace::FunctionSpaceImpl* atlas__Field__functionspace( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); return This->functionspace().get(); );
    return nullptr;
}

void atlas__Field__shapef( FieldImpl* This, int*& shape, int& rank ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); shape = const_cast<int*>( &This->shapef().front() );
                          rank                  = This->shapef().size(); );
}

void atlas__Field__host_data_int_specf( FieldImpl* This, int*& data, int& rank, int*& shapef, int*& stridesf ) {
    atlas__Field__host_data_specf( This, data, rank, shapef, stridesf );
}

void atlas__Field__host_data_long_specf( FieldImpl* This, long*& data, int& rank, int*& shapef, int*& stridesf ) {
    atlas__Field__host_data_specf( This, data, rank, shapef, stridesf );
}

void atlas__Field__host_data_float_specf( FieldImpl* This, float*& data, int& rank, int*& shapef, int*& stridesf ) {
    atlas__Field__host_data_specf( This, data, rank, shapef, stridesf );
}

void atlas__Field__host_data_double_specf( FieldImpl* This, double*& data, int& rank, int*& shapef, int*& stridesf ) {
    atlas__Field__host_data_specf( This, data, rank, shapef, stridesf );
}

void atlas__Field__device_data_int_specf( FieldImpl* This, int*& data, int& rank, int*& shapef, int*& stridesf ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); data = This->device_data<int>();
                          shapef               = const_cast<int*>( This->shapef().data() );
                          stridesf = const_cast<int*>( This->stridesf().data() ); rank = This->shapef().size(); );
}

void atlas__Field__device_data_long_specf( FieldImpl* This, long*& data, int& rank, int*& shapef, int*& stridesf ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); data = This->device_data<long>();
                          shapef               = const_cast<int*>( This->shapef().data() );
                          stridesf = const_cast<int*>( This->stridesf().data() ); rank = This->shapef().size(); );
}

void atlas__Field__device_data_float_specf( FieldImpl* This, float*& data, int& rank, int*& shapef, int*& stridesf ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); data = This->device_data<float>();
                          shapef               = const_cast<int*>( This->shapef().data() );
                          stridesf = const_cast<int*>( This->stridesf().data() ); rank = This->shapef().size(); );
}

void atlas__Field__device_data_double_specf( FieldImpl* This, double*& data, int& rank, int*& shapef, int*& stridesf ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); data = This->device_data<double>();
                          shapef               = const_cast<int*>( This->shapef().data() );
                          stridesf = const_cast<int*>( This->stridesf().data() ); rank = This->shapef().size(); );
}

int atlas__Field__host_needs_update( const FieldImpl* This ) {
    return This->hostNeedsUpdate();
}

int atlas__Field__device_needs_update( const FieldImpl* This ) {
    return This->deviceNeedsUpdate();
}

void atlas__Field__rename( FieldImpl* This, const char* name ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); This->rename( std::string( name ) ); );
}

void atlas__Field__set_levels( FieldImpl* This, int levels ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); This->set_levels( levels ); );
}

void atlas__Field__set_functionspace( FieldImpl* This, const functionspace::FunctionSpaceImpl* functionspace ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); ASSERT( functionspace ); This->set_functionspace( functionspace ); );
}

void atlas__Field__clone_to_device( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); );
    This->cloneToDevice();
}

void atlas__Field__clone_from_device( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); );
    This->cloneFromDevice();
}

void atlas__Field__sync_host_device( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); );
    This->syncHostDevice();
}

void atlas__Field__set_dirty( FieldImpl* This, int value ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); );
    This->set_dirty( value );
}

int atlas__Field__dirty( FieldImpl* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); );
    return This->dirty();
}

void atlas__Field__halo_exchange( FieldImpl* This, int on_device ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); );
    return This->haloExchange( on_device );
}
}

// ------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
