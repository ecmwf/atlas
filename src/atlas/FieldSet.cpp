// (C) Copyright 1996-2014 ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation nor
// does it submit to any jurisdiction.


#include <iostream>
#include <stdexcept>
#include <sstream>

#include "atlas/FieldSet.hpp"
#include "atlas/Field.hpp"

namespace atlas {

FieldSet::FieldSet(const std::string& name) :
  name_(name)
{
}

void FieldSet::add_field(Field& field)
{
  index_[field.name()] = fields_.size();
  fields_.push_back( &field );
}

Field& FieldSet::field(const std::string& name)
{
  try {
    return *fields_[ index_.at(name) ]; 
  }
  catch( std::out_of_range& e ) {
    std::stringstream msg;
    msg << "Could not find field \"" << name << "\" in fieldset \"" << name_ << "\"";
    throw std::out_of_range(msg.str());
  }
}

Field& FieldSet::field(size_t idx)
{ 
  return *fields_[ idx ];
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

FieldSet* atlas__FieldSet__new (char* name) {
  return new FieldSet( std::string(name) );
}

void atlas__FieldSet__delete (FieldSet* This) {
  delete This;
}

void atlas__FieldSet__add_field (FieldSet* This, Field* field) {
  This->add_field(*field);
}

int atlas__FieldSet__size(FieldSet* This) {
  return This->size();
}

void atlas__FieldSet__fields (FieldSet* This, Field** &fields, int& nb_fields)
{
  nb_fields = This->fields().size();
  fields = &This->fields()[0];
}

// ------------------------------------------------------------------


} // namespace atlas

