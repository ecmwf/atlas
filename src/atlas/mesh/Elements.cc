/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/Elements.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/library/config.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace mesh {

//-----------------------------------------------------------------------------

void Elements::rebuild() {
    size_     = hybrid_elements_->elements_size_[type_idx_];
    nb_nodes_ = hybrid_elements_->element_type( type_idx_ ).nb_nodes();
    nb_edges_ = hybrid_elements_->element_type( type_idx_ ).nb_edges();
    begin_    = hybrid_elements_->elements_begin_[type_idx_];
    end_      = hybrid_elements_->elements_begin_[type_idx_ + 1];
}

Elements::Elements( HybridElements& elements, idx_t type_idx ) :
    owns_( false ), hybrid_elements_( &elements ), type_idx_( type_idx ) {
    rebuild();
}

Elements::Elements( ElementType* element_type, idx_t nb_elements, const std::vector<idx_t>& node_connectivity ) :
    owns_( true ) {
    hybrid_elements_ = new HybridElements();
    type_idx_        = hybrid_elements_->add( element_type, nb_elements, node_connectivity.data() );
    rebuild();
}

Elements::Elements( ElementType* element_type, idx_t nb_elements, const idx_t node_connectivity[],
                    bool fortran_array ) :
    owns_( true ) {
    hybrid_elements_ = new HybridElements();
    type_idx_        = hybrid_elements_->add( element_type, nb_elements, node_connectivity, fortran_array );
    rebuild();
}

Elements::~Elements() {
    if ( owns_ ) {
        delete hybrid_elements_;
    }
}

const std::string& Elements::name() const {
    return hybrid_elements_->element_type( type_idx_ ).name();
}

template <>
array::LocalView<const double, 1> Elements::view( const Field& field ) const {
    return array::make_host_view<double, 1>( field ).slice( array::Range{begin(), begin() + size()} );
}

template <>
array::LocalView<const float, 1> Elements::view( const Field& field ) const {
    return array::make_host_view<float, 1>( field ).slice( array::Range{begin(), begin() + size()} );
}

template <>
array::LocalView<const int, 1> Elements::view( const Field& field ) const {
    return array::make_host_view<int, 1>( field ).slice( array::Range{begin(), begin() + size()} );
}

template <>
array::LocalView<const long, 1> Elements::view( const Field& field ) const {
    return array::make_host_view<long, 1>( field ).slice( array::Range{begin(), begin() + size()} );
}

template <>
array::LocalView<const double, 2> Elements::view( const Field& field ) const {
    return array::make_host_view<double, 2>( field ).slice( array::Range{begin(), begin() + size()},
                                                            array::Range::all() );
}

template <>
array::LocalView<const float, 2> Elements::view( const Field& field ) const {
    return array::make_host_view<float, 2>( field ).slice( array::Range{begin(), begin() + size()},
                                                           array::Range::all() );
}

template <>
array::LocalView<const int, 2> Elements::view( const Field& field ) const {
    return array::make_host_view<int, 2>( field ).slice( array::Range{begin(), begin() + size()}, array::Range::all() );
}

template <>
array::LocalView<const long, 2> Elements::view( const Field& field ) const {
    return array::make_host_view<long, 2>( field ).slice( array::Range{begin(), begin() + size()},
                                                          array::Range::all() );
}

// ----------------------------------------------------------------------------

template <>
array::LocalView<double, 1> Elements::view( Field& field ) const {
    return array::make_host_view<double, 1>( field ).slice( array::Range{begin(), begin() + size()} );
}

template <>
array::LocalView<float, 1> Elements::view( Field& field ) const {
    return array::make_host_view<float, 1>( field ).slice( array::Range{begin(), begin() + size()} );
}

template <>
array::LocalView<int, 1> Elements::view( Field& field ) const {
    return array::make_host_view<int, 1>( field ).slice( array::Range{begin(), begin() + size()} );
}

template <>
array::LocalView<long, 1> Elements::view( Field& field ) const {
    return array::make_host_view<long, 1>( field ).slice( array::Range{begin(), begin() + size()} );
}

template <>
array::LocalView<double, 2> Elements::view( Field& field ) const {
    return array::make_host_view<double, 2>( field ).slice( array::Range{begin(), begin() + size()},
                                                            array::Range::all() );
}

template <>
array::LocalView<float, 2> Elements::view( Field& field ) const {
    return array::make_host_view<float, 2>( field ).slice( array::Range{begin(), begin() + size()},
                                                           array::Range::all() );
}

template <>
array::LocalView<int, 2> Elements::view( Field& field ) const {
    return array::make_host_view<int, 2>( field ).slice( array::Range{begin(), begin() + size()}, array::Range::all() );
}

template <>
array::LocalView<long, 2> Elements::view( Field& field ) const {
    return array::make_host_view<long, 2>( field ).slice( array::Range{begin(), begin() + size()},
                                                          array::Range::all() );
}

idx_t Elements::add( const idx_t nb_elements ) {
    idx_t position = size();
    hybrid_elements_->insert( type_idx_, end(), nb_elements );
    return position;
}


template <>
array::LocalIndexView<idx_t, 1> Elements::indexview( Field& field ) const {
    auto local_view = array::make_host_view<idx_t, 1>( field ).slice( array::Range{begin(), begin() + size()} );
    return array::LocalIndexView<idx_t, 1>( local_view.data(), local_view.shape(), local_view.strides() );
}

template <>
array::LocalIndexView<const idx_t, 1> Elements::indexview( const Field& field ) const {
    auto local_view = array::make_host_view<idx_t, 1>( field ).slice( array::Range{begin(), begin() + size()} );
    return array::LocalIndexView<const idx_t, 1>( local_view.data(), local_view.shape(), local_view.strides() );
}

//-----------------------------------------------------------------------------

extern "C" {

void atlas__mesh__Elements__delete( Elements* This ) {
    delete This;
}

idx_t atlas__mesh__Elements__size( const Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->size();
}

idx_t atlas__mesh__Elements__begin( const Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->begin();
}

idx_t atlas__mesh__Elements__end( const Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->end();
}

BlockConnectivity* atlas__mesh__Elements__node_connectivity( Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return &This->node_connectivity();
}

BlockConnectivity* atlas__mesh__Elements__edge_connectivity( Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return &This->edge_connectivity();
}

BlockConnectivity* atlas__mesh__Elements__cell_connectivity( Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return &This->cell_connectivity();
}

int atlas__mesh__Elements__has_field( const Elements* This, char* name ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->has_field( std::string( name ) );
}

int atlas__mesh__Elements__nb_fields( const Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->nb_fields();
}

field::FieldImpl* atlas__mesh__Elements__field_by_idx( Elements* This, idx_t idx ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->field( idx ).get();
}

field::FieldImpl* atlas__mesh__Elements__field_by_name( Elements* This, char* name ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->field( std::string( name ) ).get();
}

field::FieldImpl* atlas__mesh__Elements__global_index( Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->global_index().get();
}

field::FieldImpl* atlas__mesh__Elements__remote_index( Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->remote_index().get();
}

field::FieldImpl* atlas__mesh__Elements__partition( Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->partition().get();
}

field::FieldImpl* atlas__mesh__Elements__halo( Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return This->halo().get();
}

const ElementType* atlas__mesh__Elements__element_type( const Elements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    return &This->element_type();
}

void atlas__mesh__Elements__add( Elements* This, idx_t nb_elements ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_Elements" );
    This->add( nb_elements );
}
}

//-----------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas

#undef FORTRAN_BASE
#undef TO_FORTRAN
