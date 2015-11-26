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
#include "atlas/functionspace/Spectral.h"
#include "atlas/FieldSet.h"
#include "atlas/runtime/ErrorHandling.h"

#ifdef ATLAS_HAVE_TRANS
#include "atlas/trans/Trans.h"
#endif

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------

Spectral::Spectral(const size_t truncation)
  : next::FunctionSpace(),
    truncation_(truncation),
    trans_(0)
{
}

Spectral::Spectral(trans::Trans& trans)
  : next::FunctionSpace(),
#ifdef ATLAS_HAVE_TRANS
    truncation_(trans.nsmax()),
    trans_(&trans)
#else
    truncation_(0),
    trans_(0)
#endif
{
}

Spectral::~Spectral()
{
}

size_t Spectral::nb_spectral_coefficients() const {
#ifdef ATLAS_HAVE_TRANS
  if( trans_ ) return trans_->nspec2();
#endif
  return (truncation_+1)*(truncation_+2);
}

size_t Spectral::nb_spectral_coefficients_global() const {
#ifdef ATLAS_HAVE_TRANS
  if( trans_ ) return trans_->nspec2g();
#endif
  return (truncation_+1)*(truncation_+2);
}

template <>
Field* Spectral::createField<double>(const std::string& name) const {
  Field* field = Field::create<double>(name, make_shape(nb_spectral_coefficients()) );
  field->set_functionspace(this);
  return field;
}

template <>
Field* Spectral::createField<double>(const std::string& name, size_t levels) const {
  Field* field = Field::create<double>(name, make_shape(nb_spectral_coefficients(),levels) );
  field->set_functionspace(this);
  field->set_levels(levels);
  return field;
}

template <>
Field* Spectral::createGlobalField<double>(const std::string& name) const {
  Field* field = Field::create<double>(name, make_shape(nb_spectral_coefficients_global()) );
  field->set_functionspace(this);
  return field;
}

template <>
Field* Spectral::createGlobalField<double>(const std::string& name, size_t levels) const {
  Field* field = Field::create<double>(name, make_shape(nb_spectral_coefficients_global(),levels) );
  field->set_functionspace(this);
  field->set_levels(levels);
  return field;
}

void Spectral::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const
{
  ASSERT( trans_ );
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& loc = local_fieldset[f];
    Field& glb = global_fieldset[f];
    if( loc.datatype() != DataType::str<double>() )
    {
      std::stringstream err;
      err << "Cannot gather spectral field " << loc.name() << " of datatype " << loc.datatype().str() << ".";
      err << "Only " << DataType::str<double>() << " supported.";
      throw eckit::BadValue(err.str());
    }

#ifdef ATLAS_HAVE_TRANS
    ASSERT( loc.shape(0) == nb_spectral_coefficients() );
    ASSERT( glb.shape(0) == nb_spectral_coefficients_global() );
    std::vector<int> nto(1,1);
    if( loc.rank() > 1 ) {
      nto.resize(loc.stride(0));
      for( size_t i=0; i<nto.size(); ++i ) nto[i] = 1;
    }
    trans_->gathspec(nto.size(),nto.data(),loc.data<double>(),glb.data<double>());
#else

    throw eckit::Exception("Cannot gather spectral fields because Atlas has not been compiled with TRANS support.");
#endif
  }
}
void Spectral::gather( const Field& local, Field& global ) const
{
  FieldSet local_fields;
  FieldSet global_fields;
  local_fields.add(local);
  global_fields.add(global);
  gather(local_fields,global_fields);
}

void Spectral::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const
{
  ASSERT( trans_ );
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& glb = global_fieldset[f];
    Field& loc = local_fieldset[f];
    if( loc.datatype() != DataType::str<double>() )
    {
      std::stringstream err;
      err << "Cannot scatter spectral field " << glb.name() << " of datatype " << glb.datatype().str() << ".";
      err << "Only " << DataType::str<double>() << " supported.";
      throw eckit::BadValue(err.str());
    }

#ifdef ATLAS_HAVE_TRANS
    ASSERT( loc.shape(0) == nb_spectral_coefficients() );
    ASSERT( glb.shape(0) == nb_spectral_coefficients_global() );
    std::vector<int> nfrom(1,1);
    if( loc.rank() > 1 ) {
      nfrom.resize(loc.stride(0));
      for( size_t i=0; i<nfrom.size(); ++i ) nfrom[i] = 1;
    }
    trans_->distspec(nfrom.size(),nfrom.data(),glb.data<double>(),loc.data<double>());
#else
    throw eckit::Exception("Cannot scatter spectral fields because Atlas has not been compiled with TRANS support.");
#endif

  }
}
void Spectral::scatter( const Field& global, Field& local ) const
{
  FieldSet global_fields;
  FieldSet local_fields;
  global_fields.add(global);
  local_fields.add(local);
  scatter(global_fields,local_fields);
}

std::string Spectral::checksum( const FieldSet& fieldset ) const {
  eckit::MD5 md5;
  NOTIMP;
  return md5;
}
std::string Spectral::checksum( const Field& field ) const {
  FieldSet fieldset;
  fieldset.add(field);
  return checksum(fieldset);
}

// ----------------------------------------------------------------------


extern "C"
{
Spectral* atlas__SpectralFunctionSpace__new__truncation (int truncation)
{
  ATLAS_ERROR_HANDLING(
    return new Spectral(truncation);
  );
  return 0;
}

Spectral* atlas__SpectralFunctionSpace__new__trans (trans::Trans* trans)
{
  ATLAS_ERROR_HANDLING(
    return new Spectral(*trans);
  );
  return 0;
}

void atlas__SpectralFunctionSpace__delete (Spectral* This)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    delete This;
  );
}

Field* atlas__SpectralFunctionSpace__create_field (const Spectral* This, const char* name)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField<double>(std::string(name));
  );
  return 0;
}

Field* atlas__SpectralFunctionSpace__create_field_lev (const Spectral* This, const char* name, int levels)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createField<double>(std::string(name),levels);
  );
  return 0;
}

Field* atlas__SpectralFunctionSpace__create_global_field (const Spectral* This, const char* name)
{
  ATLAS_ERROR_HANDLING (
    ASSERT(This);
    return This->createGlobalField<double>(std::string(name));
  );
  return 0;
}

Field* atlas__SpectralFunctionSpace__create_global_field_lev (const Spectral* This, const char* name, int levels)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    return This->createGlobalField<double>(std::string(name),levels);
  );
  return 0;
}

void atlas__SpectralFunctionSpace__gather (const Spectral* This, const Field* local, Field* global)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This);
    ASSERT(global);
    ASSERT(local);
    This->gather(*local,*global);
  );
}

void atlas__SpectralFunctionSpace__scatter (const Spectral* This, const Field* global, Field* local)
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

