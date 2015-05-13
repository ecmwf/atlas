/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "atlas/atlas_config.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/io/DataHandle.h"
#include "eckit/io/Buffer.h"
#include "eckit/log/Log.h"
#include "eckit/memory/ScopedPtr.h"
#include "eckit/utils/Translator.h"
#include "eckit/os/BackTrace.h"

#include "atlas/ErrorHandling.h"
#include "atlas/Field.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Grid.h"
#include "atlas/Mesh.h"
#include "atlas/Parameters.h"
#include "atlas/grids/Unstructured.h"
#include "atlas/io/PointCloud.h"
#include "atlas/util/ArrayView.h"

#include "atlas/FieldSet.h"


using namespace eckit;


namespace atlas {


//------------------------------------------------------------------------------------------------------


FieldSet::FieldSet(const std::string &name) :
  name_(name.length()? name : "untitled")
{}


void FieldSet::add_field(Field::Ptr field)
{
  index_[field->name()] = fields_.size();

  fields_.push_back( field );

  //FIXME there is a memory corruption when activating next line
//  gridset_.push_back( field->grid().self() );
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


FieldSet::FieldSet(const Field::Vector& fields) :
  fields_(fields)
{
  for( Field::Vector::const_iterator itr = fields.begin(); itr != fields.end(); ++itr )
  {
      gridset_.push_back( (*itr)->grid().self() );
  }
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

  Grid::uid_t uid = fields_[0]->grid().unique_id();

  for( size_t i = 1; i < fields_.size(); ++i )
  {
    if( fields_[i]->grid().unique_id() != uid )
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

