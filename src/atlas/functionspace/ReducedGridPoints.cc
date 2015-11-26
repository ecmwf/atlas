/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/utils/MD5.h"
#include "atlas/Mesh.h"
#include "atlas/functionspace/ReducedGridPoints.h"
#include "atlas/FieldSet.h"
#include "atlas/runtime/ErrorHandling.h"

#ifdef ATLAS_HAVE_TRANS
#include "atlas/trans/Trans.h"
#else
#error TRANS needed for functionspace::ReducedGrid
#endif

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------

ReducedGridPoints::ReducedGridPoints(const Grid& grid)
  : next::FunctionSpace(),
    grid_(grid)
{
  trans_ = new trans::Trans(grid_);
}

ReducedGridPoints::~ReducedGridPoints()
{
  delete trans_;
}

template <>
Field* ReducedGridPoints::createField<double>(const std::string& name) const {
  Field* field = Field::create<double>(name, make_shape(trans_->ngptot()) );
  field->set_functionspace(this);
  return field;
}

template <>
Field* ReducedGridPoints::createField<double>(const std::string& name, size_t levels) const {
  Field* field = Field::create<double>(name, make_shape(trans_->ngptot(),levels) );
  field->set_functionspace(this);
  field->set_levels(levels);
  return field;
}

template <>
Field* ReducedGridPoints::createGlobalField<double>(const std::string& name) const {
  Field* field = Field::create<double>(name, make_shape(trans_->ngptotg()) );
  field->set_functionspace(this);
  return field;
}

template <>
Field* ReducedGridPoints::createGlobalField<double>(const std::string& name, size_t levels) const {
  Field* field = Field::create<double>(name, make_shape(trans_->ngptotg(),levels) );
  field->set_functionspace(this);
  field->set_levels(levels);
  return field;
}

void ReducedGridPoints::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& loc = local_fieldset[f];
    Field& glb = global_fieldset[f];
    if( loc.datatype() != DataType::str<double>() )
    {
      std::stringstream err;
      err << "Cannot gather ReducedGrid field " << loc.name() << " of datatype " << loc.datatype().str() << ".";
      err << "Only " << DataType::str<double>() << " supported.";
      throw eckit::BadValue(err.str());
    }

    std::vector<int> nto(1,1);
    if( loc.rank() > 1 ) {
      nto.resize(loc.stride(0));
      for( size_t i=0; i<nto.size(); ++i ) nto[i] = 1;
    }
    trans_->gathgrid(nto.size(),nto.data(),loc.data<double>(),glb.data<double>());
  }
}
void ReducedGridPoints::gather( const Field& local, Field& global ) const
{
  FieldSet local_fields;
  FieldSet global_fields;
  local_fields.add(local);
  global_fields.add(global);
  gather(local_fields,global_fields);
}

void ReducedGridPoints::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const
{
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& glb = global_fieldset[f];
    Field& loc = local_fieldset[f];
    if( loc.datatype() != DataType::str<double>() )
    {
      std::stringstream err;
      err << "Cannot scatter ReducedGrid field " << glb.name() << " of datatype " << glb.datatype().str() << ".";
      err << "Only " << DataType::str<double>() << " supported.";
      throw eckit::BadValue(err.str());
    }

    std::vector<int> nfrom(1,1);
    if( loc.rank() > 1 ) {
      nfrom.resize(loc.stride(0));
      for( size_t i=0; i<nfrom.size(); ++i ) nfrom[i] = 1;
    }
    trans_->distgrid(nfrom.size(),nfrom.data(),glb.data<double>(),loc.data<double>());
  }
}
void ReducedGridPoints::scatter( const Field& global, Field& local ) const
{
  FieldSet global_fields;
  FieldSet local_fields;
  global_fields.add(global);
  local_fields.add(local);
  scatter(global_fields,local_fields);
}

std::string ReducedGridPoints::checksum( const FieldSet& fieldset ) const {
  eckit::MD5 md5;
  NOTIMP;
  return md5;
}
std::string ReducedGridPoints::checksum( const Field& field ) const {
  FieldSet fieldset;
  fieldset.add(field);
  return checksum(fieldset);
}

// ----------------------------------------------------------------------


extern "C"
{

ReducedGridPoints* atlas__ReducedGridFunctionSpace__new__trans (const Grid* grid)
{
  ATLAS_ERROR_HANDLING(
    return new ReducedGridPoints(*grid);
  );
  return 0;
}

void atlas__ReducedGridFunctionSpace__delete (ReducedGridPoints* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

Field* atlas__ReducedGridFunctionSpace__create_field (const ReducedGridPoints* This, const char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField<double>(std::string(name));
  );
}

Field* atlas__ReducedGridFunctionSpace__create_field_lev (const ReducedGridPoints* This, const char* name, int levels)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField<double>(std::string(name),levels);
  );
  return 0;
}

Field* atlas__ReducedGridFunctionSpace__create_global_field (const ReducedGridPoints* This, const char* name)
{
  ATLAS_ERROR_HANDLING (
    ASSERT(This);
    return This->createGlobalField<double>(std::string(name));
  );
  return 0;
}

Field* atlas__ReducedGridFunctionSpace__create_global_field_lev (const ReducedGridPoints* This, const char* name, int levels)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createGlobalField<double>(std::string(name),levels);
  );
  return 0;
}

void atlas__ReducedGridFunctionSpace__gather (const ReducedGridPoints* This, const Field* local, Field* global)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    This->gather(*local,*global);
  );
}

void atlas__ReducedGridFunctionSpace__scatter (const ReducedGridPoints* This, const Field* global, Field* local)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    This->scatter(*global,*local);
  );
}

}

} // namespace functionspace
} // namespace atlas

