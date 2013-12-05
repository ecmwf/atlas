#include "FieldSet.hpp"
#include "Field.hpp"
#include <iostream>
#include <stdexcept>
#include <sstream>
namespace ecmwf {

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

FieldSet* ecmwf__FieldSet__new (char* name) {
  return new FieldSet( std::string(name) );
}

void ecmwf__FieldSet__delete (FieldSet* This) {
  delete This;
}

void ecmwf__FieldSet__add_field (FieldSet* This, Field* field) {
  This->add_field(*field);
}

int ecmwf__FieldSet__size(FieldSet* This) {
  return This->size();
}

void ecmwf__FieldSet__fields (FieldSet* This, Field** &fields, int& nb_fields)
{
  nb_fields = This->fields().size();
  fields = &This->fields()[0];
}

// ------------------------------------------------------------------


} // namespace ecmwf

