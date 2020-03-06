/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>

#include "eckit/log/Bytes.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/library/config.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

#if ATLAS_HAVE_FORTRAN
#define FORTRAN_BASE 1
#define TO_FORTRAN +1
#else
#define FORTRAN_BASE 0
#endif

using atlas::array::ArrayView;
using atlas::array::IndexView;
using atlas::array::make_datatype;
using atlas::array::make_indexview;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

namespace {

static void set_uninitialized_fields_to_zero( HybridElements& elems, idx_t begin ) {
    ArrayView<gidx_t, 1> global_index = make_view<gidx_t, 1>( elems.global_index() );
    IndexView<idx_t, 1> remote_index  = make_indexview<idx_t, 1>( elems.remote_index() );
    ArrayView<int, 1> partition       = make_view<int, 1>( elems.partition() );
    ArrayView<int, 1> halo            = make_view<int, 1>( elems.halo() );
    ArrayView<int, 1> flags           = make_view<int, 1>( elems.flags() );

    for ( idx_t j = begin; j < elems.size(); ++j ) {
        global_index( j ) = 0;
        remote_index( j ) = 0;
        partition( j )    = 0;
        halo( j )         = 0;
        flags( j )        = 0;
    }
}
}  // namespace

//------------------------------------------------------------------------------------------------------

HybridElements::HybridElements() : size_( 0 ), elements_size_(), elements_begin_( 1, 0ul ), type_idx_() {
    add( Field( "glb_idx", make_datatype<gidx_t>(), make_shape( size() ) ) );
    add( Field( "remote_idx", make_datatype<idx_t>(), make_shape( size() ) ) );
    add( Field( "partition", make_datatype<int>(), make_shape( size() ) ) );
    add( Field( "halo", make_datatype<int>(), make_shape( size() ) ) );
    add( Field( "flags", make_datatype<int>(), make_shape( size() ) ) );
    set_uninitialized_fields_to_zero( *this, 0 );

    node_connectivity_ = &add( new Connectivity( "node" ) );
    edge_connectivity_ = &add( new Connectivity( "edge" ) );
    cell_connectivity_ = &add( new Connectivity( "cell" ) );
}

HybridElements::~HybridElements() = default;

Field HybridElements::add( const Field& field ) {
    ATLAS_ASSERT( field );
    ATLAS_ASSERT( !field.name().empty() );

    if ( has_field( field.name() ) ) {
        std::stringstream msg;
        msg << "Trying to add field '" << field.name()
            << "' to HybridElements, but HybridElements already has a field with "
               "this name.";
        throw_Exception( msg.str(), Here() );
    }
    fields_[field.name()] = field;
    return field;
}

void HybridElements::resize( idx_t size ) {
    idx_t old_size = size_;
    size_          = size;
    for ( FieldMap::iterator it = fields_.begin(); it != fields_.end(); ++it ) {
        Field& field            = it->second;
        array::ArrayShape shape = field.shape();
        shape[0]                = size_;
        field.resize( shape );
    }

    set_uninitialized_fields_to_zero( *this, old_size );
}

void HybridElements::remove_field( const std::string& name ) {
    if ( !has_field( name ) ) {
        std::stringstream msg;
        msg << "Trying to remove field `" << name
            << "' in HybridElements, but no field with this name is present in HybridElements.";
        throw_Exception( msg.str(), Here() );
    }
    fields_.erase( name );
}

const Field& HybridElements::field( const std::string& name ) const {
    if ( !has_field( name ) ) {
        std::stringstream msg;
        msg << "Trying to access field `" << name
            << "' in HybridElements, but no field with this name is present in HybridElements.";
        throw_Exception( msg.str(), Here() );
    }
    return fields_.find( name )->second;
}

Field& HybridElements::field( const std::string& name ) {
    return const_cast<Field&>( static_cast<const HybridElements*>( this )->field( name ) );
}

const Field& HybridElements::field( idx_t idx ) const {
    ATLAS_ASSERT( idx < nb_fields() );
    idx_t c( 0 );
    for ( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it ) {
        if ( idx == c ) {
            const Field& field = it->second;
            return field;
        }
        c++;
    }
    throw_Exception( "Should not be here!", Here() );
}

Field& HybridElements::field( idx_t idx ) {
    return const_cast<Field&>( static_cast<const HybridElements*>( this )->field( idx ) );
}

HybridElements::Connectivity& HybridElements::add( Connectivity* connectivity ) {
    connectivities_[connectivity->name()] = connectivity;
    return *connectivity;
}

idx_t HybridElements::add( const ElementType* element_type, idx_t nb_elements,
                           const std::vector<idx_t>& connectivity ) {
    return add( element_type, nb_elements, connectivity.data() );
}

idx_t HybridElements::add( const ElementType* element_type, idx_t nb_elements, const idx_t connectivity[] ) {
    return add( element_type, nb_elements, connectivity, false );
}

idx_t HybridElements::add( const ElementType* element_type, idx_t nb_elements, const idx_t connectivity[],
                           bool fortran_array ) {
    util::ObjectHandle<const ElementType> etype( element_type );

    idx_t old_size = size();
    idx_t new_size = old_size + nb_elements;

    idx_t nb_nodes = etype->nb_nodes();

    type_idx_.resize( new_size );

    for ( idx_t e = old_size; e < new_size; ++e ) {
        type_idx_[e] = element_types_.size();
    }

    elements_begin_.push_back( new_size );
    elements_size_.push_back( nb_elements );

    element_types_.push_back( etype );
    elements_.resize( element_types_.size() );
    for ( idx_t t = 0; t < nb_types(); ++t ) {
        if ( elements_[t] ) {
            elements_[t]->rebuild();
        }
        else {
            elements_[t] = util::ObjectHandle<Elements>( new Elements( *this, t ) );
        }
    }

    node_connectivity_->add( nb_elements, nb_nodes, connectivity, fortran_array );

    resize( new_size );
    return element_types_.size() - 1;
}

idx_t HybridElements::add( const ElementType* element_type, idx_t nb_elements ) {
    util::ObjectHandle<const ElementType> etype( element_type );

    idx_t old_size = size();
    idx_t new_size = old_size + nb_elements;

    idx_t nb_nodes = etype->nb_nodes();

    type_idx_.resize( new_size );

    for ( idx_t e = old_size; e < new_size; ++e ) {
        type_idx_[e] = element_types_.size();
    }

    elements_begin_.emplace_back( new_size );
    elements_size_.emplace_back( nb_elements );

    element_types_.emplace_back( etype );
    elements_.emplace_back( new Elements( *this, elements_.size() ) );

    node_connectivity_->add( nb_elements, nb_nodes );
    resize( new_size );
    return element_types_.size() - 1;
}

idx_t HybridElements::add( const Elements& elems ) {
    bool fortran_array = true;
    return add( &elems.element_type(), elems.size(), elems.node_connectivity().data(), fortran_array );
}

const std::string& HybridElements::name( idx_t elem_idx ) const {
    return element_types_[type_idx_[elem_idx]]->name();
}

idx_t HybridElements::elemtype_nb_nodes( idx_t elem_idx ) const {
    return element_type( type_idx( elem_idx ) ).nb_nodes();
}

idx_t HybridElements::elemtype_nb_edges( idx_t elem_idx ) const {
    return element_type( type_idx( elem_idx ) ).nb_edges();
}

void HybridElements::insert( idx_t type_idx, idx_t position, idx_t nb_elements ) {
    type_idx_.insert( type_idx_.begin() + position, nb_elements, type_idx );
    elements_size_[type_idx] += nb_elements;
    for ( idx_t jtype = type_idx + 1; jtype < nb_types() + 1; ++jtype ) {
        elements_begin_[jtype] += nb_elements;
    }
    for ( idx_t t = 0; t < nb_types(); ++t ) {
        elements_[t]->rebuild();
    }
    node_connectivity_->insert( position, nb_elements, element_types_[type_idx]->nb_nodes() );

    size_ += nb_elements;
    for ( FieldMap::iterator it = fields_.begin(); it != fields_.end(); ++it ) {
        Field& field = it->second;
        field.insert( position, nb_elements );
    }
}

//-----------------------------------------------------------------------------

void HybridElements::clear() {
    resize( 0 );
    for ( ConnectivityMap::iterator it = connectivities_.begin(); it != connectivities_.end(); ++it ) {
        it->second->clear();
    }
    size_ = 0;
    elements_size_.clear();
    elements_begin_.resize( 1 );
    elements_begin_[0] = 0;
    element_types_.clear();
    type_idx_.clear();
    elements_.clear();
}

//-----------------------------------------------------------------------------

void HybridElements::updateDevice() const {
    std::for_each( fields_.begin(), fields_.end(), []( const FieldMap::value_type& v ) { v.second.updateDevice(); } );
}

void HybridElements::updateHost() const {
    std::for_each( fields_.begin(), fields_.end(), []( const FieldMap::value_type& v ) { v.second.updateHost(); } );
}

void HybridElements::syncHostDevice() const {
    std::for_each( fields_.begin(), fields_.end(), []( const FieldMap::value_type& v ) { v.second.syncHostDevice(); } );
}

size_t HybridElements::footprint() const {
    size_t size = sizeof( *this );
    for ( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it ) {
        size += ( *it ).second.footprint();
    }
    for ( ConnectivityMap::const_iterator it = connectivities_.begin(); it != connectivities_.end(); ++it ) {
        size += ( *it ).second->footprint();
    }
    size += elements_size_.capacity() * sizeof( idx_t );
    size += elements_begin_.capacity() * sizeof( idx_t );

    size += metadata_.footprint();

    return size;
}

//-----------------------------------------------------------------------------

extern "C" {
HybridElements* atlas__mesh__HybridElements__create() {
    return new HybridElements();
}

void atlas__mesh__HybridElements__delete( HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    delete This;
}

MultiBlockConnectivity* atlas__mesh__HybridElements__node_connectivity( HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return &This->node_connectivity();
}

MultiBlockConnectivity* atlas__mesh__HybridElements__edge_connectivity( HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return &This->edge_connectivity();
}

MultiBlockConnectivity* atlas__mesh__HybridElements__cell_connectivity( HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return &This->cell_connectivity();
}

idx_t atlas__mesh__HybridElements__size( const HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return This->size();
}

void atlas__mesh__HybridElements__add_elements( HybridElements* This, ElementType* elementtype, idx_t nb_elements ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    This->add( elementtype, nb_elements );
}

void atlas__mesh__HybridElements__add_elements_with_nodes( HybridElements* This, ElementType* elementtype,
                                                           idx_t nb_elements, idx_t node_connectivity[],
                                                           int fortran_array ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    ATLAS_ASSERT( elementtype != nullptr, "Cannot access uninitialised atlas_mesh_ElementType" );
    This->add( elementtype, nb_elements, node_connectivity, fortran_array );
}

int atlas__mesh__HybridElements__has_field( const HybridElements* This, char* name ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return This->has_field( std::string( name ) );
}

int atlas__mesh__HybridElements__nb_fields( const HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return This->nb_fields();
}

int atlas__mesh__HybridElements__nb_types( const HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return This->nb_types();
}

field::FieldImpl* atlas__mesh__HybridElements__field_by_idx( HybridElements* This, idx_t idx ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return This->field( idx ).get();
}

field::FieldImpl* atlas__mesh__HybridElements__field_by_name( HybridElements* This, char* name ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return This->field( std::string( name ) ).get();
}

field::FieldImpl* atlas__mesh__HybridElements__global_index( HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return This->global_index().get();
}

field::FieldImpl* atlas__mesh__HybridElements__remote_index( HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return This->remote_index().get();
}

field::FieldImpl* atlas__mesh__HybridElements__partition( HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return This->partition().get();
}

field::FieldImpl* atlas__mesh__HybridElements__halo( HybridElements* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return This->halo().get();
}

Elements* atlas__mesh__HybridElements__elements( HybridElements* This, idx_t idx ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    return &This->elements( idx );
}

void atlas__mesh__HybridElements__add_field( HybridElements* This, field::FieldImpl* field ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_mesh_HybridElements" );
    This->add( field );
}
}

//-----------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas

#undef FORTRAN_BASE
