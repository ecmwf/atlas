/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "StructuredColumnsInterface.h"

#include "atlas/field/FieldSet.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/distribution/DistributionImpl.h"
#include "atlas/runtime/ErrorHandling.h"

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------------
// Fortran interfaces
// ----------------------------------------------------------------------------

namespace detail {
struct StructuredColumnsFortranAccess {
    detail::StructuredColumns::Map2to1& ij2gp_;
    StructuredColumnsFortranAccess( const detail::StructuredColumns& fs ) :
        ij2gp_( const_cast<detail::StructuredColumns&>( fs ).ij2gp_ ) {}
};
}  // namespace detail


extern "C" {

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid(
    const Grid::Implementation* grid, const eckit::Configuration* config ) {
    ATLAS_ERROR_HANDLING( return new detail::StructuredColumns( Grid( grid ), grid::Partitioner(), *config ); );
    return nullptr;
}

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_dist(
    const Grid::Implementation* grid, const grid::DistributionImpl* dist, const eckit::Configuration* config ) {
    ATLAS_ERROR_HANDLING( return new detail::StructuredColumns( Grid( grid ), grid::Distribution( dist ), *config ); );
    return nullptr;
}


void atlas__functionspace__StructuredColumns__delete( detail::StructuredColumns* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); delete This; );
}

field::FieldImpl* atlas__fs__StructuredColumns__create_field( const detail::StructuredColumns* This,
                                                              const eckit::Configuration* options ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); field::FieldImpl * field; {
        Field f = This->createField( *options );
        field   = f.get();
        field->attach();
    } field->detach();
                          return field; );
    return nullptr;
}

void atlas__functionspace__StructuredColumns__gather( const detail::StructuredColumns* This,
                                                      const field::FieldImpl* local, field::FieldImpl* global ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); ASSERT( global ); ASSERT( local ); const Field l( local ); Field g( global );
                          This->gather( l, g ); );
}

void atlas__functionspace__StructuredColumns__scatter( const detail::StructuredColumns* This,
                                                       const field::FieldImpl* global, field::FieldImpl* local ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); ASSERT( global ); ASSERT( local ); const Field g( global ); Field l( local );
                          This->scatter( g, l ); );
}

void atlas__fs__StructuredColumns__halo_exchange_field( const detail::StructuredColumns* This,
                                                        const field::FieldImpl* field ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); ASSERT( field ); Field f( field ); This->haloExchange( f ); );
}

void atlas__fs__StructuredColumns__halo_exchange_fieldset( const detail::StructuredColumns* This,
                                                           const field::FieldSetImpl* fieldset ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); ASSERT( fieldset ); FieldSet f( fieldset ); This->haloExchange( f ); );
}

void atlas__fs__StructuredColumns__checksum_fieldset( const detail::StructuredColumns* This,
                                                      const field::FieldSetImpl* fieldset, char*& checksum, idx_t& size,
                                                      int& allocated ) {
    ASSERT( This );
    ASSERT( fieldset );
    ATLAS_ERROR_HANDLING( std::string checksum_str( This->checksum( fieldset ) );
                          size = static_cast<idx_t>( checksum_str.size() ); checksum = new char[size + 1];
                          allocated = true; strcpy( checksum, checksum_str.c_str() ); );
}

void atlas__fs__StructuredColumns__checksum_field( const detail::StructuredColumns* This, const field::FieldImpl* field,
                                                   char*& checksum, idx_t& size, int& allocated ) {
    ASSERT( This );
    ASSERT( field );
    ATLAS_ERROR_HANDLING( std::string checksum_str( This->checksum( field ) );
                          size = static_cast<idx_t>( checksum_str.size() ); checksum = new char[size + 1];
                          allocated = true; strcpy( checksum, checksum_str.c_str() ); );
}

void atlas__fs__StructuredColumns__index_host( const detail::StructuredColumns* This, idx_t*& data, idx_t& i_min,
                                               idx_t& i_max, idx_t& j_min, idx_t& j_max ) {
    ASSERT( This );
    auto _This = detail::StructuredColumnsFortranAccess{*This};
    ATLAS_ERROR_HANDLING( data = _This.ij2gp_.data_.data(); i_min = _This.ij2gp_.i_min_ + 1;
                          i_max = _This.ij2gp_.i_max_ + 1; j_min = _This.ij2gp_.j_min_ + 1;
                          j_max                                  = _This.ij2gp_.j_max_ + 1; );
}

idx_t atlas__fs__StructuredColumns__j_begin( const detail::StructuredColumns* This ) {
    return This->j_begin() + 1;
}
idx_t atlas__fs__StructuredColumns__j_end( const detail::StructuredColumns* This ) {
    return This->j_end();
}
idx_t atlas__fs__StructuredColumns__i_begin( const detail::StructuredColumns* This, idx_t j ) {
    return This->i_begin( j - 1 ) + 1;
}
idx_t atlas__fs__StructuredColumns__i_end( const detail::StructuredColumns* This, idx_t j ) {
    return This->i_end( j - 1 );
}
idx_t atlas__fs__StructuredColumns__j_begin_halo( const detail::StructuredColumns* This ) {
    return This->j_begin_halo() + 1;
}
idx_t atlas__fs__StructuredColumns__j_end_halo( const detail::StructuredColumns* This ) {
    return This->j_end_halo();
}
idx_t atlas__fs__StructuredColumns__i_begin_halo( const detail::StructuredColumns* This, idx_t j ) {
    return This->i_begin_halo( j - 1 ) + 1;
}
idx_t atlas__fs__StructuredColumns__i_end_halo( const detail::StructuredColumns* This, idx_t j ) {
    return This->i_end_halo( j - 1 );
}

field::FieldImpl* atlas__fs__StructuredColumns__xy( const detail::StructuredColumns* This ) {
    return This->xy().get();
}

field::FieldImpl* atlas__fs__StructuredColumns__partition( const detail::StructuredColumns* This ) {
    return This->partition().get();
}

field::FieldImpl* atlas__fs__StructuredColumns__global_index( const detail::StructuredColumns* This ) {
    return This->global_index().get();
}

field::FieldImpl* atlas__fs__StructuredColumns__index_i( const detail::StructuredColumns* This ) {
    return This->index_i().get();
}

field::FieldImpl* atlas__fs__StructuredColumns__index_j( const detail::StructuredColumns* This ) {
    return This->index_j().get();
}
}

// ----------------------------------------------------------------------------

}  // namespace functionspace
}  // namespace atlas
