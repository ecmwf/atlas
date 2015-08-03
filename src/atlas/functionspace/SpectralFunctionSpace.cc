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

SpectralFunctionSpace::SpectralFunctionSpace(const std::string& name, const trans::Trans& trans )
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

// ----------------------------------------------------------------------

SpectralColumnFunctionSpace::SpectralColumnFunctionSpace(const std::string& name, const size_t truncation, const size_t nb_levels)
  : SpectralFunctionSpace(name,truncation),
    nb_levels_(nb_levels)
{
}

SpectralColumnFunctionSpace::SpectralColumnFunctionSpace(const std::string& name, const trans::Trans& trans, const size_t nb_levels )
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

