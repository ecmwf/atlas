/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_SpectralFunctionSpace_h
#define atlas_functionspace_SpectralFunctionSpace_h

#include "atlas/atlas_defines.h"
#include "atlas/FunctionSpace.h"

namespace atlas { namespace trans { class Trans; } }
namespace atlas { class Field; }
namespace atlas { class FieldSet; }

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class Spectral : public FunctionSpace
{
public:

  Spectral(const size_t truncation);

  Spectral(trans::Trans& );

  virtual ~Spectral();

  virtual std::string name() const { return "Spectral"; }

  /// @brief Create a spectral field
  template <typename DATATYPE> Field* createField(const std::string& name) const;
  template <typename DATATYPE> Field* createField(const std::string& name, size_t levels) const;

  /// @brief Create a global spectral field
  template <typename DATATYPE> Field* createGlobalField(const std::string& name) const;
  template <typename DATATYPE> Field* createGlobalField(const std::string& name, size_t levels) const;

  void gather( const FieldSet&, FieldSet& ) const;
  void gather( const Field&, Field& ) const;

  void scatter( const FieldSet&, FieldSet& ) const;
  void scatter( const Field&, Field& ) const;

  std::string checksum( const FieldSet& ) const;
  std::string checksum( const Field& ) const;

public: // methods

  size_t nb_spectral_coefficients() const;
  size_t nb_spectral_coefficients_global() const;

private: // data

  size_t truncation_;

  trans::Trans* trans_;

};

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
#define Trans trans::Trans
extern "C"
{
  Spectral* atlas__SpectralFunctionSpace__new__truncation (int truncation);
  Spectral* atlas__SpectralFunctionSpace__new__trans (Trans* trans);
  void atlas__SpectralFunctionSpace__delete (Spectral* This);

  Field* atlas__SpectralFunctionSpace__create_field (const Spectral* This, const char* name);

  Field* atlas__SpectralFunctionSpace__create_field_lev (const Spectral* This, const char* name, int levels);

  Field* atlas__SpectralFunctionSpace__create_global_field (const Spectral* This, const char* name);

  Field* atlas__SpectralFunctionSpace__create_global_field_lev (const Spectral* This, const char* name, int levels);

  void atlas__SpectralFunctionSpace__gather (const Spectral* This, const Field* local, Field* global);

  void atlas__SpectralFunctionSpace__scatter (const Spectral* This, const Field* global, Field* local);


}
#undef Trans

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_SpectralFunctionSpace_h
