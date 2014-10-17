/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <iostream>
#include <stdexcept>
#include <sstream>

#include "eckit/exception/Exceptions.h"
#include "atlas/FieldGroup.h"
#include "atlas/Field.h"

namespace atlas {

//------------------------------------------------------------------------------------------------------

FieldGroup::FieldGroup(const std::string& name) :
	name_(name)
{
}

void FieldGroup::add_field(Field& field)
{
	index_[field.name()] = fields_.size();
	fields_.push_back( &field );
}

Field& FieldGroup::field(const std::string& name)
{
	if( has_field(name) )
	{
		return *fields_[ index_.at(name) ];
	}
	else
	{
		std::stringstream msg;
		msg << "Could not find field \"" << name << "\" in fieldset \"" << name_ << "\"";
		throw eckit::OutOfRange(msg.str(),Here());
	}
}

Field& FieldGroup::field(size_t idx)
{
	return *fields_[ idx ];
}

//------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

FieldGroup* atlas__FieldSet__new (char* name) {
	return new FieldGroup( std::string(name) );
}

void atlas__FieldSet__delete (FieldGroup* This) {
	delete This;
}

void atlas__FieldSet__add_field (FieldGroup* This, Field* field) {
	This->add_field(*field);
}

int atlas__FieldSet__has_field (FieldGroup* This, char* name) {
	return This->has_field( std::string(name) );
}

int atlas__FieldSet__size(FieldGroup* This) {
	return This->size();
}

void atlas__FieldSet__fields (FieldGroup* This, Field** &fields, int& nb_fields)
{
	nb_fields = This->fields().size();
	fields = &This->fields()[0];
}

Field* atlas__FieldSet__field_by_name (FieldGroup* This, char* name) {
	return &This->field( std::string(name) );
}

Field* atlas__FieldSet__field_by_idx (FieldGroup* This, int idx) {
	return &This->field( idx );
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

