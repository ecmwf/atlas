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

#include <eckit/exception/Exceptions.h>
#include <eckit/filesystem/PathName.h>
#include <eckit/io/DataHandle.h>
#include <eckit/io/Buffer.h>
#include <eckit/log/Log.h>
#include <eckit/memory/ScopedPtr.h>
#include <eckit/utils/Translator.h>
#include "eckit/os/BackTrace.h"

#ifdef ECKIT_HAVE_GRIB
#include <eckit/grib/GribField.h>
#include <eckit/grib/GribFieldSet.h>
#include <eckit/grib/GribHandle.h>
#include <eckit/grib/GribParams.h>
#endif

#include "atlas/Field.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Parameters.h"
#include "atlas/util/ArrayView.h"
#include "atlas/io/PointCloud.h"

#include "atlas/Mesh.h"
#include "atlas/Grid.h"
#include "atlas/grids/Unstructured.h"

#include "atlas/FieldSet.h"


using namespace eckit;


namespace atlas {


//------------------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------


FieldSet::FieldSet(const std::string &name) :
  name_(name.length()? name : "untitled")
{}


void FieldSet::add_field(const Field& field)
{
  index_[field.name()] = fields_.size();
  fields_.push_back( Field::Ptr(const_cast< Field* >(&field)) );
}

bool FieldSet::has_field(const std::string& name) const
{
  return index_.count(name);
}

Field& FieldSet::field(const std::string& name) const
{
  if (!has_field(name))
  {
    const std::string msg("Could not find field \"" + name + "\" in fieldset \"" + name_ + "\"");
    throw eckit::OutOfRange(msg,Here());
  }
  return *fields_[ index_.at(name) ];
}


FieldSet::FieldSet(const PathName& path) :
  fields_(),
  grid_()
{
  // read start of file (using buffer) to determine file format
  Buffer buff(64);
  DataHandle* dh = path.fileHandle();
  dh->openForRead();
  size_t len = dh->read(buff,buff.size());
  ASSERT(len);
  dh->close();


  //TODO: make and use a factory


  if ((len>=4) && (0==strncmp(buff,"GRIB",4)))
  {
#ifdef ECKIT_HAVE_GRIB
    // attempt to read GRIB format
    using grib::GribField;
    using grib::GribFieldSet;
    using grib::GribHandle;

    GribFieldSet gribfs(path);
    if (!gribfs.size())
      return;

    fields_.reserve(gribfs.size());
    for (size_t fidx=0; fidx<gribfs.size(); ++fidx)
    {
      GribField* gf = const_cast< GribField* >( gribfs.get(fidx) );
      ASSERT(gf);

      GribHandle* gh = gf->getHandle();
      ASSERT(gh);

      // (grid_ built on first call)
      fields_.push_back( create_field(*gh) );

      gf->release();
    }

    ASSERT(grid_);
    ASSERT(haveSameGrid());
#else
    std::stringstream stream;
    stream << "eckit was built without GRIB support, cannot construct FieldSet from GRIB path " << path << "\n";
    stream << eckit::BackTrace::dump();
    throw eckit::Exception(stream.str(), Here());
#endif
    return;
  }


  // attempt to read PointCloud format
  if ((len>=10) && (0==strncmp(buff,"PointCloud",10)))
  {
    ScopedPtr< grids::Unstructured > grid(io::PointCloud::read(path));
    grid->mesh();
    NOTIMP;
    return;
  }
}


FieldSet::FieldSet(const Buffer& buf) :
	fields_(),
	grid_()
{
#ifdef ECKIT_HAVE_GRIB
  grib::GribHandle gh(buf);

  // (grid_ built inside)
  fields_.push_back( create_field(gh) );

  ASSERT( grid_ );
  ASSERT( haveSameGrid() );
#else
  eckit::Log::error() << eckit::BackTrace::dump() << std::endl;
  throw eckit::Exception("eckit was built without GRIB support, cannot construct FieldSet from GRIB buffer", Here());
#endif
}


FieldSet::FieldSet(const DataHandle&)
{
  NOTIMP;
}


FieldSet::FieldSet(const Grid::Ptr grid, const std::vector<std::string>& nfields) :
	fields_(),
	grid_(grid)
{
	ASSERT( grid_ );

	Mesh& mesh = grid_->mesh();

	FunctionSpace& nodes = mesh.function_space( "nodes" );
    fields_.reserve(nfields.size());
    for( size_t i = 0; i < nfields.size(); ++i )
    {
		Field& f = nodes.create_field<double>(nfields[i],1);

		ASSERT( grid->uid() == f.grid().uid() );

		fields_.push_back( Field::Ptr( &f ) );
    }

  ASSERT( haveSameGrid() );
}


FieldSet::FieldSet(const Field::Vector& fields) :
	fields_(fields),
	grid_()
{
	if( !fields_.empty() )
		grid_.reset( fields_[0]->grid().self() );

	ASSERT( grid_ );
  ASSERT( haveSameGrid() );
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


#if 0
void FieldSet::add_field(Field& field)
{
  index_[field.name()] = fields_.size();
  fields_.push_back( Field::Ptr(&field) );
}
#endif


#if 0
Field& FieldSet::field(const std::string& name)
{
  if(!has_field(name))
  {
    const std::string msg("FieldSet" + (name_.length()? " \"" + name_ + "\"" : "") + ": cannot find field \"" + name + "\"");
    throw eckit::OutOfRange(msg,Here());
  }
  return *fields_[ index_.at(name) ];
}
#endif


#if 0
const Field& FieldSet::field(const std::string& name) const
{
  if(!has_field(name))
  {
    const std::string msg("FieldSet" + (name_.length()? " \"" + name_ + "\"" : "") + ": cannot find field \"" + name + "\"");
    throw eckit::OutOfRange(msg,Here());
  }
  return *fields_[ index_.at(name) ];
}
#endif


Field::Ptr FieldSet::create_field(grib::GribHandle& gh)
{
#ifdef ECKIT_HAVE_GRIB
  if(!grid_)
  {
    grib::GribParams* gp = grib::GribParams::create(gh);
    ASSERT( gp );
    grid_.reset( Grid::create( *gp ) );
  }
  else // check grid is the same
  {
    if( gh.geographyHash() != grid_->hash() )
      throw eckit::UserError("GRIB fields don't match grid within FieldSet", Here() );
  }

  Mesh& mesh = grid_->mesh();
    FunctionSpace&  nodes  = mesh.function_space( "nodes" );

    // get name for this field
  std::string sname = gh.shortName() + "_" + Translator<size_t,std::string>()( fields_.size() );

    // get values

    size_t nvalues = gh.getDataValuesSize();

    // create the field

    if( nodes.shape(0) != nvalues )
        throw SeriousBug( "Size of field in GRIB does not match Grid", Here() );

  Field& f = nodes.create_field<double>(sname,1);

  gh.getDataValues( f.data<double>(), nvalues );

  f.grib( gh.clone() );

  return f.self();

#else
  throw eckit::Exception("eckit was built without GRIB support\n  --> Cannot create field from GribHandle", Here());
  return Field::Ptr();
#endif
}


bool FieldSet::haveSameGrid() const
{
  if( fields_.empty() ) return true;

  bool result = true;

  std::string uid = grid_->uid();

  for( size_t i = 0; i < fields_.size(); ++i )
  {
    if( fields_[i]->grid().uid() != uid )
    {
      Log::error() << "fields_["<< i <<"] (" << fields_[i]->name() << ") doesn't match the same grid as FieldSet @ " << Here() << std::endl;
      result = false;
    }
  }

  return result;
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
  FieldSet* fset = new FieldSet( std::string(name) );
  fset->name() = name;
  return fset;
}

void atlas__FieldSet__delete(FieldSet* This)
{ delete This; }


void atlas__FieldSet__fields (FieldSet* This, Field** &fields, int &nb_fields)
{
  nb_fields = This->fields().size();

  if (fields!=NULL)
    throw eckit::SeriousBug("provided return pointer is not NULL (memory leak)");

  fields = nb_fields ? __private_get_raw_fields_ptr(This).data() : NULL;
}


void   atlas__FieldSet__add_field     (FieldSet* This, Field* field) 
{ ASSERT(This != NULL); This->add_field(*field); }

int    atlas__FieldSet__has_field     (FieldSet* This, char* name)   
{ ASSERT(This != NULL); return This->has_field( std::string(name) ); }

int    atlas__FieldSet__size          (FieldSet* This)               
{ ASSERT(This != NULL); return This->size(); }

Field* atlas__FieldSet__field_by_name (FieldSet* This, char* name)   
{ ASSERT(This != NULL); return &This->field( std::string(name) ); }

Field* atlas__FieldSet__field_by_idx  (FieldSet* This, int idx)      
{ ASSERT(This != NULL); return &This->operator[](idx); }

}
//-----------------------------------------------------------------------------


}  // namespace atlas

