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

namespace atlas {
namespace functionspace {

// ----------------------------------------------------------------------

SpectralFunctionSpace::SpectralFunctionSpace(const std::string& name, const size_t truncation)
  : next::FunctionSpace(name),
    truncation_(truncation)
{
}

SpectralFunctionSpace::~SpectralFunctionSpace() {}

size_t SpectralFunctionSpace::nspec2g() const {
  return (truncation_+1)*(truncation_+2);
}

template<> Field* SpectralFunctionSpace::create_field<float>(const std::string& name) {
  return new field::FieldT<float>(name, make_shape(nspec2g()) );
}

template<> Field* SpectralFunctionSpace::create_field<double>(const std::string& name) {
  return new field::FieldT<double>(name, make_shape(nspec2g()) );
}

// ----------------------------------------------------------------------

SpectralColumnFunctionSpace::SpectralColumnFunctionSpace(const std::string& name, const size_t truncation, const size_t levels)
  : next::FunctionSpace(name),
    truncation_(truncation),
    levels_(levels)
{
}

SpectralColumnFunctionSpace::~SpectralColumnFunctionSpace() {}

size_t SpectralColumnFunctionSpace::nspec2g() const {
  return (truncation_+1)*(truncation_+2);
}

template<> Field* SpectralColumnFunctionSpace::create_field<float>(const std::string& name) {
  return new field::FieldT<float>(name, make_shape(nspec2g(),levels_) );
}

template<> Field* SpectralColumnFunctionSpace::create_field<double>(const std::string& name) {
  return new field::FieldT<double>(name, make_shape(nspec2g(),levels_) );
}


} // namespace functionspace
} // namespace atlas

