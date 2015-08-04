/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/Mesh.h"
#include "atlas/functionspace/SpectralFunctionSpace.h"
#include "atlas/field/FieldT.h"
#include "atlas/FieldSet.h"

#ifdef ATLAS_HAVE_TRANS
#include "atlas/trans/Trans.h"
#endif

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------

SpectralFunctionSpace::SpectralFunctionSpace(const std::string& name, const size_t truncation)
  : next::FunctionSpace(name),
    truncation_(truncation),
    trans_(0)
{
}

SpectralFunctionSpace::SpectralFunctionSpace(const std::string& name, trans::Trans& trans)
  : next::FunctionSpace(name),
#ifdef ATLAS_HAVE_TRANS
    truncation_(trans.nsmax()),
    trans_(&trans)
#else
    truncation_(0),
    trans_(0)
#endif
{
}

SpectralFunctionSpace::~SpectralFunctionSpace()
{
}

size_t SpectralFunctionSpace::nb_spectral_coefficients() const {
#ifdef ATLAS_HAVE_TRANS
  if( trans_ ) return trans_->nspec2();
#endif
  return (truncation_+1)*(truncation_+2);
}

size_t SpectralFunctionSpace::nb_spectral_coefficients_global() const {
#ifdef ATLAS_HAVE_TRANS
  if( trans_ ) return trans_->nspec2g();
#endif
  return (truncation_+1)*(truncation_+2);
}

Field* SpectralFunctionSpace::createField(const std::string& name) {
  return new field::FieldT<double>(name, make_shape(nb_spectral_coefficients()) );
}

Field* SpectralFunctionSpace::createGlobalField(const std::string& name) {
  return new field::FieldT<double>(name, make_shape(nb_spectral_coefficients_global()) );
}

void SpectralFunctionSpace::gather( const FieldSet& local_fieldset, FieldSet& global_fieldset ) const
{
  ASSERT( trans_ );
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& loc = local_fieldset[f];
    Field& glb = global_fieldset[f];
    if( loc.datatype() != DataType::datatype<double>() )
    {
      std::stringstream err;
      err << "Cannot gather spectral field " << loc.name() << " of datatype " << loc.datatype() << ".";
      err << "Only " << DataType::datatype<double>() << " supported.";
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
void SpectralFunctionSpace::gather( const Field& local, Field& global ) const
{
  FieldSet local_fields;
  FieldSet global_fields;
  local_fields.add(local);
  global_fields.add(global);
  gather(local_fields,global_fields);
}

void SpectralFunctionSpace::scatter( const FieldSet& global_fieldset, FieldSet& local_fieldset ) const
{
  ASSERT( trans_ );
  ASSERT(local_fieldset.size() == global_fieldset.size());

  for( size_t f=0; f<local_fieldset.size(); ++f ) {

    const Field& glb = global_fieldset[f];
    Field& loc = local_fieldset[f];
    if( loc.datatype() != DataType::datatype<double>() )
    {
      std::stringstream err;
      err << "Cannot scatter spectral field " << glb.name() << " of datatype " << glb.datatype() << ".";
      err << "Only " << DataType::datatype<double>() << " supported.";
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
void SpectralFunctionSpace::scatter( const Field& global, Field& local ) const
{
  FieldSet global_fields;
  FieldSet local_fields;
  global_fields.add(global);
  local_fields.add(local);
  scatter(global_fields,local_fields);
}

std::string SpectralFunctionSpace::checksum( const FieldSet& fieldset ) const {
  eckit::MD5 md5;
  NOTIMP;
  return md5;
}
std::string SpectralFunctionSpace::checksum( const Field& field ) const {
  FieldSet fieldset;
  fieldset.add(field);
  return checksum(fieldset);
}


// ----------------------------------------------------------------------

SpectralColumnFunctionSpace::SpectralColumnFunctionSpace(const std::string& name, const size_t truncation, const size_t nb_levels)
  : SpectralFunctionSpace(name,truncation),
    nb_levels_(nb_levels)
{
}

SpectralColumnFunctionSpace::SpectralColumnFunctionSpace(const std::string& name, trans::Trans& trans, const size_t nb_levels)
  : SpectralFunctionSpace(name,trans),
    nb_levels_(nb_levels)
{
}

SpectralColumnFunctionSpace::~SpectralColumnFunctionSpace() {}

size_t SpectralColumnFunctionSpace::nb_levels() const {
  return nb_levels_;
}

Field* SpectralColumnFunctionSpace::createField(const std::string& name) {
  return new field::FieldT<double>(name, make_shape(nb_spectral_coefficients(),nb_levels_) );
}

Field* SpectralColumnFunctionSpace::createGlobalField(const std::string& name) {
  return new field::FieldT<double>(name, make_shape(nb_spectral_coefficients_global(),nb_levels_) );
}

} // namespace functionspace
} // namespace atlas

