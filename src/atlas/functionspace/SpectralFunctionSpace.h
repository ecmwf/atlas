/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_SpectralFunctionSpace_h
#define atlas_functionspace_SpectralFunctionSpace_h

#include "atlas/FunctionSpace.h"

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class SpectralFunctionSpace : public next::FunctionSpace
{
public:

  SpectralFunctionSpace(const std::string& name, const size_t truncation);

  virtual ~SpectralFunctionSpace();

  /// @brief Create a spectral field
  template< typename DATATYPE >
  Field* create_field(const std::string& name);

private: // methods

  size_t nspec2g() const;

private: // data

  size_t truncation_;

};

// -------------------------------------------------------------------

class SpectralColumnFunctionSpace : public next::FunctionSpace
{
public:

  SpectralColumnFunctionSpace(const std::string& name, const size_t truncation, const size_t levels);

  virtual ~SpectralColumnFunctionSpace();

  /// @brief Create a spectral field
  template< typename DATATYPE >
  Field* create_field(const std::string& name);

private: // methods

  size_t nspec2g() const;

private: // data

  size_t truncation_;
  size_t levels_;

};


} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_SpectralFunctionSpace_h
