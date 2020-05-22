/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cstring>

#include "StructuredColumnsInterface.h"

#include "atlas/field/FieldSet.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/distribution/DistributionImpl.h"
#include "atlas/runtime/Exception.h"

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
    return new detail::StructuredColumns( Grid( grid ), grid::Partitioner(), *config );
}

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_dist(
    const Grid::Implementation* grid, const grid::DistributionImpl* dist, const eckit::Configuration* config ) {
    return new detail::StructuredColumns( Grid( grid ), grid::Distribution( dist ), *config );
}

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_dist_vert(
    const Grid::Implementation* grid, const grid::DistributionImpl* dist, const Vertical* vert,
    const eckit::Configuration* config ) {
    return new detail::StructuredColumns( Grid( grid ), grid::Distribution( dist ), *vert, *config );
}
const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_dist_config(
    const Grid::Implementation* grid, const grid::DistributionImpl* dist, const eckit::Configuration* config ) {
    return new detail::StructuredColumns( Grid( grid ), grid::Distribution( dist ), *config );
}

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_part(
    const Grid::Implementation* grid, const PartitionerImpl* partitioner, const eckit::Configuration* config ) {
    return new detail::StructuredColumns( Grid( grid ), grid::Partitioner( partitioner ), *config );
}

const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid_part_vert(
    const Grid::Implementation* grid, const PartitionerImpl* partitioner, const Vertical* vert,
    const eckit::Configuration* config ) {
    return new detail::StructuredColumns( Grid( grid ), *vert, grid::Partitioner( partitioner ), *config );
}

void atlas__functionspace__StructuredColumns__gather_field( const detail::StructuredColumns* This,
                                                            const field::FieldImpl* local, field::FieldImpl* global ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_functionspace_StructuredColumns" );
    ATLAS_ASSERT( global != nullptr, "Cannot access uninitialised atlas_Field" );
    ATLAS_ASSERT( local != nullptr, "Cannot access uninitialised atlas_Field" );
    const Field l( local );
    Field g( global );
    This->gather( l, g );
}

void atlas__functionspace__StructuredColumns__scatter_field( const detail::StructuredColumns* This,
                                                             const field::FieldImpl* global, field::FieldImpl* local ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_functionspace_StructuredColumns" );
    ATLAS_ASSERT( global != nullptr, "Cannot access uninitialised atlas_Field" );
    ATLAS_ASSERT( local != nullptr, "Cannot access uninitialised atlas_Field" );
    const Field g( global );
    Field l( local );
    This->scatter( g, l );
}

void atlas__functionspace__StructuredColumns__gather_fieldset( const detail::StructuredColumns* This,
                                                               const field::FieldSetImpl* local,
                                                               field::FieldSetImpl* global ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_functionspace_StructuredColumns" );
    ATLAS_ASSERT( global != nullptr, "Cannot access uninitialised atlas_FieldSet" );
    ATLAS_ASSERT( local != nullptr, "Cannot access uninitialised atlas_FieldSet" );
    const FieldSet l( local );
    FieldSet g( global );
    This->gather( l, g );
}

void atlas__functionspace__StructuredColumns__scatter_fieldset( const detail::StructuredColumns* This,
                                                                const field::FieldSetImpl* global,
                                                                field::FieldSetImpl* local ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_functionspace_StructuredColumns" );
    ATLAS_ASSERT( global != nullptr, "Cannot access uninitialised atlas_FieldSet" );
    ATLAS_ASSERT( local != nullptr, "Cannot access uninitialised atlas_FieldSet" );
    const FieldSet g( global );
    FieldSet l( local );
    This->scatter( g, l );
}

void atlas__fs__StructuredColumns__checksum_fieldset( const detail::StructuredColumns* This,
                                                      const field::FieldSetImpl* fieldset, char*& checksum, idx_t& size,
                                                      int& allocated ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_functionspace_StructuredColumns" );
    ATLAS_ASSERT( fieldset != nullptr, "Cannot access uninitialised atlas_FieldSet" );
    std::string checksum_str( This->checksum( fieldset ) );
    size      = static_cast<idx_t>( checksum_str.size() );
    checksum  = new char[size + 1];
    allocated = true;
    strcpy( checksum, checksum_str.c_str() );
}

void atlas__fs__StructuredColumns__checksum_field( const detail::StructuredColumns* This, const field::FieldImpl* field,
                                                   char*& checksum, idx_t& size, int& allocated ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_functionspace_StructuredColumns" );
    ATLAS_ASSERT( field != nullptr, "Cannot access uninitialised atlas_Field" );
    std::string checksum_str( This->checksum( field ) );
    size      = static_cast<idx_t>( checksum_str.size() );
    checksum  = new char[size + 1];
    allocated = true;
    strcpy( checksum, checksum_str.c_str() );
}

void atlas__fs__StructuredColumns__index_host( const detail::StructuredColumns* This, idx_t*& data, idx_t& i_min,
                                               idx_t& i_max, idx_t& j_min, idx_t& j_max ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_functionspace_StructuredColumns" );
    auto _This = detail::StructuredColumnsFortranAccess{*This};
    data       = _This.ij2gp_.data_.data();
    i_min      = _This.ij2gp_.i_min_ + 1;
    i_max      = _This.ij2gp_.i_max_ + 1;
    j_min      = _This.ij2gp_.j_min_ + 1;
    j_max      = _This.ij2gp_.j_max_ + 1;
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

field::FieldImpl* atlas__fs__StructuredColumns__z( const detail::StructuredColumns* This ) {
    ATLAS_ASSERT( This != nullptr );
    field::FieldImpl* field;
    {
        Field f = This->z();
        field   = f.get();
        field->attach();
    }
    field->detach();
    return field;
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

idx_t atlas__fs__StructuredColumns__size( const detail::StructuredColumns* This ) {
    return This->size();
}

idx_t atlas__fs__StructuredColumns__sizeOwned( const detail::StructuredColumns* This ) {
    return This->sizeOwned();
}

idx_t atlas__fs__StructuredColumns__levels( const detail::StructuredColumns* This ) {
    return This->levels();
}

const GridImpl* atlas__fs__StructuredColumns__grid( const detail::StructuredColumns* This ) {
    return This->grid().get();
}
}


// ----------------------------------------------------------------------------

}  // namespace functionspace
}  // namespace atlas
