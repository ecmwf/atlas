/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/Field.h"
#include "atlas/Grid.h"
#include "atlas/FieldSet.h"

namespace atlas {

//------------------------------------------------------------------------------------------------------


FieldSet::FieldSet(const std::string &name) :
  name_(name.length()? name : "untitled")
{}

void FieldSet::add(Field& field)
{
  index_[field.name()] = fields_.size();
  fields_.push_back( field.self() );
}

void FieldSet::add_field(Field::Ptr field)
{
  index_[field->name()] = fields_.size();
  fields_.push_back( field );
}


bool FieldSet::has_field(const std::string& name) const
{
  return index_.count(name);
}


Field& FieldSet::field(const std::string& name) const
{
  if (!has_field(name))
  {
    const std::string msg("FieldSet" + (name_.length()? " \"" + name_ + "\"" : "") + ": cannot find field \"" + name + "\"");
    throw eckit::OutOfRange(msg,Here());
  }
  return *fields_[ index_.at(name) ];
}


FieldSet::FieldSet(const std::vector< Field::Ptr >& fields) :
  fields_(fields)
{
}


std::vector< std::string > FieldSet::field_names() const
{
  std::vector< std::string > ret;
  if (fields_.size())
    ret.reserve(fields_.size());

  for (size_t i=0; i<fields_.size(); ++i)
    ret.push_back(fields_[i]->name());

  return ret;
}


bool FieldSet::have_same_grid() const
{
  if( fields_.empty() )
    return true;

  Grid::uid_t uid = fields_[0]->grid().uniqueId();

  for( size_t i = 1; i < fields_.size(); ++i )
  {
    if( fields_[i]->grid().uniqueId() != uid )
      return false;
  }

  return true;
}


std::vector<Field*>& __private_get_raw_fields_ptr (FieldSet* This)
{
  This->fields_raw_ptr_.resize( This->size() );
  for( int f=0; f<This->size(); ++f )
    This->fields_raw_ptr_[f] = This->fields()[f].get();
  return This->fields_raw_ptr_;
}

//-----------------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"{


FieldSet* atlas__FieldSet__new (char* name)
{
  ATLAS_ERROR_HANDLING(
    FieldSet* fset = new FieldSet( std::string(name) );
    fset->name() = name;
    return fset;
  );
  return NULL;
}


void atlas__FieldSet__delete(FieldSet* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT( This!= NULL );
    delete This;
  );
}


void atlas__FieldSet__fields (FieldSet* This, Field** &fields, int &nb_fields)
{
  ATLAS_ERROR_HANDLING(
    nb_fields = This->fields().size();

    if (fields!=NULL)
      throw eckit::SeriousBug("provided return pointer is not NULL (memory leak)");

    fields = nb_fields ? __private_get_raw_fields_ptr(This).data() : NULL;
  );
}


void   atlas__FieldSet__add_field     (FieldSet* This, Field* field)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This != NULL);
    This->add_field(field->self());
  );
}

int    atlas__FieldSet__has_field     (FieldSet* This, char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This != NULL);
    return This->has_field( std::string(name) );
  );
  return 0;
}

int    atlas__FieldSet__size          (FieldSet* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This != NULL);
    return This->size();
  );
  return 0;
}

Field* atlas__FieldSet__field_by_name (FieldSet* This, char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This != NULL);
    return &This->field( std::string(name) );
  );
  return NULL;
}

Field* atlas__FieldSet__field_by_idx  (FieldSet* This, int idx)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This != NULL);
    return &This->operator[](idx);
  );
  return NULL;
}


}
//-----------------------------------------------------------------------------


}  // namespace atlas

