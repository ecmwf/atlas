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
#include "atlas/Grid.h"
#include "atlas/Mesh.h"
#include "atlas/functionspace/ReducedGridPoints.h"
#include "atlas/FieldSet.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/util/Checksum.h"

#ifdef ATLAS_HAVE_TRANS
#include "atlas/trans/Trans.h"
#endif

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------

ReducedGridPoints::ReducedGridPoints(const Grid& grid)
  : next::FunctionSpace(),
    grid_(grid)
{
#ifdef ATLAS_HAVE_TRANS
  trans_ = new trans::Trans(grid_);
#endif
}

ReducedGridPoints::~ReducedGridPoints()
{
#ifdef ATLAS_HAVE_TRANS
  delete trans_;
#endif
}

template <>
Field* ReducedGridPoints::createField<double>(const std::string& name) const {
#ifdef ATLAS_HAVE_TRANS
  Field* field = Field::create<double>(name, make_shape(trans_->ngptot()) );
  field->set_functionspace(this);
  return field;
#else
  eckit::NotImplemented("ReducedGridPoints::createField currently relies on ATLAS_HAVE_TRANS",Here());
  Field* field = Field::create<double>(name, make_shape(grid_.npts()) );
  field->set_functionspace(this);
  return field;
#endif
}

template <>
Field* ReducedGridPoints::createField<double>(const std::string& name, size_t levels) const {
#ifdef ATLAS_HAVE_TRANS
  Field* field = Field::create<double>(name, make_shape(trans_->ngptot(),levels) );
  field->set_functionspace(this);
  field->set_levels(levels);
  return field;
#else
  eckit::NotImplemented("ReducedGridPoints::createField currently relies on ATLAS_HAVE_TRANS",Here());
  Field* field = Field::create<double>(name, make_shape(grid_.npts(),levels) );
  field->set_functionspace(this);
  return field;
#endif
}

template <>
Field* ReducedGridPoints::createGlobalField<double>(const std::string& name) const {
#ifdef ATLAS_HAVE_TRANS
  Field* field = Field::create<double>(name, make_shape(trans_->ngptotg()) );
  field->set_functionspace(this);
  return field;
#else
  eckit::NotImplemented("ReducedGridPoints::createGlobalField currently relies on ATLAS_HAVE_TRANS",Here());
  Field* field = Field::create<double>(name, make_shape(grid_.npts()) );
  field->set_functionspace(this);
  return field;
#endif
}

template <>
Field* ReducedGridPoints::createGlobalField<double>(const std::string& name, size_t levels) const {
#ifdef ATLAS_HAVE_TRANS
  Field* field = Field::create<double>(name, make_shape(trans_->ngptotg(),levels) );
  field->set_functionspace(this);
  field->set_levels(levels);
  return field;
#else
  eckit::NotImplemented("ReducedGridPoints::createGlobalField currently relies on ATLAS_HAVE_TRANS",Here());
  Field* field = Field::create<double>(name, make_shape(grid_.npts(),levels) );
  field->set_functionspace(this);
  return field;
#endif
}

void ReducedGridPoints::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const
{
#ifdef ATLAS_HAVE_TRANS
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
#else
  eckit::NotImplemented("ReducedGridPoints::gather currently relies on ATLAS_HAVE_TRANS",Here());
#endif
}
void ReducedGridPoints::gather( const Field& local, Field& global ) const
{
#ifdef ATLAS_HAVE_TRANS
  FieldSet local_fields;
  FieldSet global_fields;
  local_fields.add(local);
  global_fields.add(global);
  gather(local_fields,global_fields);
#else
  eckit::NotImplemented("ReducedGridPoints::gather currently relies on ATLAS_HAVE_TRANS",Here());
#endif
}

void ReducedGridPoints::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const
{
#ifdef ATLAS_HAVE_TRANS
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
#else
  eckit::NotImplemented("ReducedGridPoints::scatter currently relies on ATLAS_HAVE_TRANS",Here());
#endif
}

void ReducedGridPoints::scatter( const Field& global, Field& local ) const
{
#ifdef ATLAS_HAVE_TRANS
  FieldSet global_fields;
  FieldSet local_fields;
  global_fields.add(global);
  local_fields.add(local);
  scatter(global_fields,local_fields);
#else
  eckit::NotImplemented("ReducedGridPoints::scatter currently relies on ATLAS_HAVE_TRANS",Here());
#endif
}

std::string ReducedGridPoints::checksum( const FieldSet& fieldset ) const {
  eckit::MD5 md5;
  NOTIMP;
  return md5;
}
std::string ReducedGridPoints::checksum( const Field& field ) const {
  // FieldSet fieldset;
  // fieldset.add(field);  
  // return checksum(fieldset);
  eckit::Log::warning() << "Only local checksum implemented" << std::endl;
  std::stringstream resultss;
  resultss << atlas::checksum(field.data<double>(),field.size());
  return resultss.str();
}

// ----------------------------------------------------------------------


extern "C"
{

ReducedGridPoints* atlas__functionspace__ReducedGridPoints__new__grid (const Grid* grid)
{
  ATLAS_ERROR_HANDLING(
    return new ReducedGridPoints(*grid);
  );
  return 0;
}

void atlas__functionspace__ReducedGridPoints__delete (ReducedGridPoints* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

Field* atlas__functionspace__ReducedGridPoints__create_field (const ReducedGridPoints* This, const char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField<double>(std::string(name));
  );
  return 0;
}

Field* atlas__functionspace__ReducedGridPoints__create_field_lev (const ReducedGridPoints* This, const char* name, int levels)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField<double>(std::string(name),levels);
  );
  return 0;
}

Field* atlas__functionspace__ReducedGridPoints__create_gfield (const ReducedGridPoints* This, const char* name)
{
  ATLAS_ERROR_HANDLING (
    ASSERT(This);
    return This->createGlobalField<double>(std::string(name));
  );
  return 0;
}

Field* atlas__functionspace__ReducedGridPoints__create_gfield_lev (const ReducedGridPoints* This, const char* name, int levels)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createGlobalField<double>(std::string(name),levels);
  );
  return 0;
}

void atlas__functionspace__ReducedGridPoints__gather (const ReducedGridPoints* This, const Field* local, Field* global)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    This->gather(*local,*global);
  );
}

void atlas__functionspace__ReducedGridPoints__scatter (const ReducedGridPoints* This, const Field* global, Field* local)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    This->scatter(*global,*local);
  );
}

void atlas__fs__ReducedGridPoints__checksum_fieldset(const ReducedGridPoints* This, const FieldSet* fieldset, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(fieldset);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(*fieldset));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

void atlas__fs__ReducedGridPoints__checksum_field(const ReducedGridPoints* This, const Field* field, char* &checksum, int &size, int &allocated)
{
  ASSERT(This);
  ASSERT(field);
  ATLAS_ERROR_HANDLING(
    std::string checksum_str (This->checksum(*field));
    size = checksum_str.size();
    checksum = new char[size+1]; allocated = true;
    strcpy(checksum,checksum_str.c_str());
  );
}

}

} // namespace functionspace
} // namespace atlas

