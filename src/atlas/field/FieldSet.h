/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file FieldSet.h
/// @author Willem Deconinck
/// @author Pedro Maciel
/// @date Jan 2015

#pragma once

#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "atlas/field/Field.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"

namespace atlas {

class FieldSet;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace field {

/**
 * @brief Represents a set of fields, where order is preserved
 */
class FieldSetImpl : public util::Object {
public:  // types
    using iterator       = std::vector<Field>::iterator;
    using const_iterator = std::vector<Field>::const_iterator;

public:  // methods
    /// Constructs an empty FieldSet
    FieldSetImpl( const std::string& name = "untitled" );

    idx_t size() const { return static_cast<idx_t>( fields_.size() ); }
    bool empty() const { return !fields_.size(); }

    void clear();

    const std::string& name() const { return name_; }
    std::string& name() { return name_; }

    const Field& operator[]( const idx_t& i ) const { return field( i ); }
    Field& operator[]( const idx_t& i ) { return field( i ); }

    const Field& operator[]( const std::string& name ) const { return field( name ); }
    Field& operator[]( const std::string& name ) { return field( name ); }

    const Field& field( const idx_t& i ) const {
        if ( i >= size() )
            throw_OutOfRange( "fieldset", i, size(), Here() );
        return fields_[i];
    }
    Field& field( const idx_t& i ) {
        if ( i >= size() )
            throw_OutOfRange( "fieldset", i, size(), Here() );
        return fields_[i];
    }

    std::vector<std::string> field_names() const;

    Field add( const Field& );

    bool has_field( const std::string& name ) const;

    Field& field( const std::string& name ) const;

    iterator begin() { return fields_.begin(); }
    iterator end() { return fields_.end(); }
    const_iterator begin() const { return fields_.begin(); }
    const_iterator end() const { return fields_.end(); }
    const_iterator cbegin() const { return fields_.begin(); }
    const_iterator cend() const { return fields_.end(); }

    void haloExchange( bool on_device = false ) const;
    void adjointHaloExchange( bool on_device = false ) const;
    void set_dirty( bool = true ) const;

protected:                                // data
    std::vector<Field> fields_;           ///< field storage
    std::string name_;                    ///< internal name
    std::map<std::string, idx_t> index_;  ///< name-to-index map, to refer fields by name
};

class FieldImpl;

// C wrapper interfaces to C++ routines
extern "C" {
FieldSetImpl* atlas__FieldSet__new( char* name );
void atlas__FieldSet__delete( FieldSetImpl* This );
void atlas__FieldSet__add_field( FieldSetImpl* This, FieldImpl* field );
int atlas__FieldSet__has_field( const FieldSetImpl* This, char* name );
idx_t atlas__FieldSet__size( const FieldSetImpl* This );
FieldImpl* atlas__FieldSet__field_by_name( FieldSetImpl* This, char* name );
FieldImpl* atlas__FieldSet__field_by_idx( FieldSetImpl* This, idx_t idx );
void atlas__FieldSet__set_dirty( FieldSetImpl* This, int value );
void atlas__FieldSet__halo_exchange( FieldSetImpl* This, int on_device );
}

}  // namespace field
#endif

//---------------------------------------------------------------------------------------------------------------------

/**
 * @brief Represents a set of fields, where order is preserved
 */
class FieldSet : DOXYGEN_HIDE( public util::ObjectHandle<field::FieldSetImpl> ) {
public:  // types
    using iterator       = Implementation::iterator;
    using const_iterator = Implementation::const_iterator;

public:  // methods
    using Handle::Handle;
    FieldSet();
    FieldSet( const std::string& name );
    FieldSet( const Field& );

    idx_t size() const { return get()->size(); }
    bool empty() const { return get()->empty(); }

    void clear() { get()->clear(); }

    const std::string& name() const { return get()->name(); }
    std::string& name() { return get()->name(); }

    const Field& operator[]( const idx_t& i ) const { return get()->operator[]( i ); }
    Field& operator[]( const idx_t& i ) { return get()->operator[]( i ); }

    const Field& operator[]( const std::string& name ) const { return get()->operator[]( name ); }
    Field& operator[]( const std::string& name ) { return get()->operator[]( name ); }

    const Field& operator[]( const char* name ) const { return get()->operator[]( name ); }
    Field& operator[]( const char* name ) { return get()->operator[]( name ); }

    const Field& field( const idx_t& i ) const { return get()->field( i ); }
    Field& field( const idx_t& i ) { return get()->field( i ); }

    std::vector<std::string> field_names() const { return get()->field_names(); }

    Field add( const Field& field ) { return get()->add( field ); }

    bool has_field( const std::string& name ) const { return get()->has_field( name ); }

    Field& field( const std::string& name ) const { return get()->field( name ); }

    iterator begin() { return get()->begin(); }
    iterator end() { return get()->end(); }
    const_iterator begin() const { return get()->begin(); }
    const_iterator end() const { return get()->end(); }
    const_iterator cbegin() const { return get()->begin(); }
    const_iterator cend() const { return get()->end(); }

    void haloExchange( bool on_device = false ) const { get()->haloExchange( on_device ); }
    void set_dirty( bool = true ) const;
};

}  // namespace atlas
