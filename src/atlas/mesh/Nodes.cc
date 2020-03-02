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

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"

using atlas::array::make_datatype;
using atlas::array::make_shape;

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

Nodes::Nodes() : size_( 0 ) {
    global_index_ = add( Field( "glb_idx", make_datatype<gidx_t>(), make_shape( size() ) ) );
    remote_index_ = add( Field( "remote_idx", make_datatype<idx_t>(), make_shape( size() ) ) );
    partition_    = add( Field( "partition", make_datatype<int>(), make_shape( size() ) ) );
    xy_           = add( Field( "xy", make_datatype<double>(), make_shape( size(), 2 ) ) );
    xy_.set_variables( 2 );
    lonlat_ = add( Field( "lonlat", make_datatype<double>(), make_shape( size(), 2 ) ) );
    lonlat_.set_variables( 2 );
    ghost_ = add( Field( "ghost", make_datatype<int>(), make_shape( size() ) ) );
    flags_ = add( Field( "flags", make_datatype<int>(), make_shape( size() ) ) );
    halo_  = add( Field( "halo", make_datatype<int>(), make_shape( size() ) ) );

    edge_connectivity_ = &add( new Connectivity( "edge" ) );
    cell_connectivity_ = &add( new Connectivity( "cell" ) );
}

Nodes::Connectivity& Nodes::add( Connectivity* connectivity ) {
    connectivities_[connectivity->name()] = connectivity;
    return *connectivity;
}

Field Nodes::add( const Field& field ) {
    ATLAS_ASSERT( field );
    ATLAS_ASSERT( !field.name().empty() );

    if ( has_field( field.name() ) ) {
        std::stringstream msg;
        msg << "Trying to add field '" << field.name() << "' to Nodes, but Nodes already has a field with this name.";
        throw_Exception( msg.str(), Here() );
    }
    fields_[field.name()] = field;
    return field;
}

void Nodes::remove_field( const std::string& name ) {
    if ( !has_field( name ) ) {
        std::stringstream msg;
        msg << "Trying to remove field `" << name << "' in Nodes, but no field with this name is present in Nodes.";
        throw_Exception( msg.str(), Here() );
    }
    fields_.erase( name );
}

const Field& Nodes::field( const std::string& name ) const {
    if ( !has_field( name ) ) {
        std::stringstream msg;
        msg << "Trying to access field `" << name << "' in Nodes, but no field with this name is present in Nodes.";
        throw_Exception( msg.str(), Here() );
    }
    return fields_.find( name )->second;
}

Field& Nodes::field( const std::string& name ) {
    return const_cast<Field&>( static_cast<const Nodes*>( this )->field( name ) );
}

void Nodes::resize( idx_t size ) {
    if ( size != size_ ) {
        idx_t previous_size = size_;
        size_               = size;
        for ( FieldMap::iterator it = fields_.begin(); it != fields_.end(); ++it ) {
            Field& field            = it->second;
            array::ArrayShape shape = field.shape();
            shape[0]                = size_;
            field.resize( shape );
        }

        auto glb_idx = array::make_view<gidx_t, 1>( global_index() );
        auto part    = array::make_view<int, 1>( partition() );
        auto flag    = array::make_view<int, 1>( flags() );
        auto _halo   = array::make_view<int, 1>( halo() );

        const int mpi_rank = static_cast<int>( mpi::rank() );
        for ( idx_t n = previous_size; n < size_; ++n ) {
            glb_idx( n ) = 1 + n;
            part( n )    = mpi_rank;
            flag( n )    = 0;
            _halo( n )   = std::numeric_limits<int>::max();
        }
    }
}

const Field& Nodes::field( idx_t idx ) const {
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
    static Field f;
    return f;
}
Field& Nodes::field( idx_t idx ) {
    return const_cast<Field&>( static_cast<const Nodes*>( this )->field( idx ) );
}

void Nodes::print( std::ostream& os ) const {
    os << "Nodes[\n";
    os << "\t size=" << size() << ",\n";
    os << "\t fields=\n";
    for ( idx_t i = 0; i < nb_fields(); ++i ) {
        os << "\t\t" << field( i );
        if ( i != nb_fields() - 1 ) {
            os << ",";
        }
        os << "\n";
    }
    os << "]";
}

size_t Nodes::footprint() const {
    size_t size = sizeof( *this );
    for ( FieldMap::const_iterator it = fields_.begin(); it != fields_.end(); ++it ) {
        size += ( *it ).second.footprint();
    }
    for ( ConnectivityMap::const_iterator it = connectivities_.begin(); it != connectivities_.end(); ++it ) {
        size += ( *it ).second->footprint();
    }
    size += metadata_.footprint();
    return size;
}

const IrregularConnectivity& Nodes::connectivity( const std::string& name ) const {
    if ( ( connectivities_.find( name ) == connectivities_.end() ) ) {
        std::stringstream msg;
        msg << "Trying to access connectivity `" << name
            << "' in Nodes, but no connectivity with this name is present in "
               "Nodes.";
        throw_Exception( msg.str(), Here() );
    }
    return *connectivities_.find( name )->second;
}
IrregularConnectivity& Nodes::connectivity( const std::string& name ) {
    if ( ( connectivities_.find( name ) == connectivities_.end() ) ) {
        std::stringstream msg;
        msg << "Trying to access connectivity `" << name
            << "' in Nodes, but no connectivity with this name is present in "
               "Nodes.";
        throw_Exception( msg.str(), Here() );
    }
    return *connectivities_.find( name )->second;
}

void Nodes::updateDevice() const {
    std::for_each( fields_.begin(), fields_.end(), []( const FieldMap::value_type& v ) { v.second.updateDevice(); } );
}

void Nodes::updateHost() const {
    std::for_each( fields_.begin(), fields_.end(), []( const FieldMap::value_type& v ) { v.second.updateHost(); } );
}

void Nodes::syncHostDevice() const {
    std::for_each( fields_.begin(), fields_.end(), []( const FieldMap::value_type& v ) { v.second.syncHostDevice(); } );
}

//-----------------------------------------------------------------------------

extern "C" {

Nodes* atlas__mesh__Nodes__create() {
    return new Nodes();
}

void atlas__mesh__Nodes__delete( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
}

idx_t atlas__mesh__Nodes__size( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return This->size();
}
void atlas__mesh__Nodes__resize( Nodes* This, idx_t size ) {
    ATLAS_ASSERT( This != nullptr );
    This->resize( size );
}
idx_t atlas__mesh__Nodes__nb_fields( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return This->nb_fields();
}

void atlas__mesh__Nodes__add_field( Nodes* This, field::FieldImpl* field ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( field != nullptr );
    This->add( field );
}

void atlas__mesh__Nodes__remove_field( Nodes* This, char* name ) {
    ATLAS_ASSERT( This != nullptr );
    This->remove_field( std::string( name ) );
}

int atlas__mesh__Nodes__has_field( Nodes* This, char* name ) {
    ATLAS_ASSERT( This != nullptr );
    return This->has_field( std::string( name ) );
}

field::FieldImpl* atlas__mesh__Nodes__field_by_name( Nodes* This, char* name ) {
    ATLAS_ASSERT( This != nullptr );
    return This->field( std::string( name ) ).get();
}

field::FieldImpl* atlas__mesh__Nodes__field_by_idx( Nodes* This, idx_t idx ) {
    ATLAS_ASSERT( This != nullptr );
    return This->field( idx ).get();
}

util::Metadata* atlas__mesh__Nodes__metadata( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return &This->metadata();
}

void atlas__mesh__Nodes__str( Nodes* This, char*& str, int& size ) {
    ATLAS_ASSERT( This != nullptr );
    std::stringstream ss;
    ss << *This;
    std::string s = ss.str();
    size          = static_cast<int>( s.size() );
    str           = new char[size + 1];
    strcpy( str, s.c_str() );
}

IrregularConnectivity* atlas__mesh__Nodes__edge_connectivity( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return &This->edge_connectivity();
}

IrregularConnectivity* atlas__mesh__Nodes__cell_connectivity( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return &This->cell_connectivity();
}

IrregularConnectivity* atlas__mesh__Nodes__connectivity( Nodes* This, char* name ) {
    ATLAS_ASSERT( This != nullptr );
    return &This->connectivity( std::string( name ) );
}

void atlas__mesh__Nodes__add_connectivity( Nodes* This, IrregularConnectivity* connectivity ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( connectivity != nullptr );
    This->add( connectivity );
}

field::FieldImpl* atlas__mesh__Nodes__xy( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return This->xy().get();
}

field::FieldImpl* atlas__mesh__Nodes__lonlat( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return This->lonlat().get();
}

field::FieldImpl* atlas__mesh__Nodes__global_index( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return This->global_index().get();
}

field::FieldImpl* atlas__mesh__Nodes__remote_index( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return This->remote_index().get();
}

field::FieldImpl* atlas__mesh__Nodes__partition( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return This->partition().get();
}

field::FieldImpl* atlas__mesh__Nodes__ghost( Nodes* This ) {
    ATLAS_ASSERT( This != nullptr );
    return This->ghost().get();
}
}

}  // namespace mesh
}  // namespace atlas
