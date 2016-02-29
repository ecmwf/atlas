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
#include "atlas/functionspace/FunctionSpace.h"

namespace atlas {
namespace field {
    class Field;
    class FieldSet;
} }

namespace atlas {
namespace numerics {
namespace trans {
    class Trans;
} } }

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class Spectral : public FunctionSpace
{
public:

  Spectral(const size_t truncation);

  Spectral(numerics::trans::Trans& );

  virtual ~Spectral();

  virtual std::string name() const { return "Spectral"; }

  /// @brief Create a spectral field
  template <typename DATATYPE> field::Field* createField(
    const std::string& name) const;
  template <typename DATATYPE> field::Field* createField(
    const std::string& name,
    size_t levels) const;

  /// @brief Create a global spectral field
  template <typename DATATYPE> field::Field* createGlobalField(
    const std::string& name) const;
  template <typename DATATYPE> field::Field* createGlobalField(
    const std::string& name, size_t levels) const;

  void gather( const field::FieldSet&, field::FieldSet& ) const;
  void gather( const field::Field&, field::Field& ) const;

  void scatter( const field::FieldSet&, field::FieldSet& ) const;
  void scatter( const field::Field&, field::Field& ) const;

  std::string checksum( const field::FieldSet& ) const;
  std::string checksum( const field::Field& ) const;

public: // methods

  size_t nb_spectral_coefficients() const;
  size_t nb_spectral_coefficients_global() const;

private: // data

  size_t truncation_;

  numerics::trans::Trans* trans_;

};

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
#define Trans numerics::trans::Trans
#define field_Field field::Field
extern "C"
{
  Spectral* atlas__SpectralFunctionSpace__new__truncation (int truncation);
  Spectral* atlas__SpectralFunctionSpace__new__trans (Trans* trans);
  void atlas__SpectralFunctionSpace__delete (Spectral* This);
  field_Field* atlas__SpectralFunctionSpace__create_field(const Spectral* This, const char* name);
  field_Field* atlas__SpectralFunctionSpace__create_field_lev(const Spectral* This, const char* name, int levels);
  field_Field* atlas__SpectralFunctionSpace__create_global_field(const Spectral* This, const char* name);
  field_Field* atlas__SpectralFunctionSpace__create_global_field_lev(const Spectral* This, const char* name, int levels);
  void atlas__SpectralFunctionSpace__gather(const Spectral* This, const field_Field* local, field_Field* global);
  void atlas__SpectralFunctionSpace__scatter(const Spectral* This, const field_Field* global, field_Field* local);
}
#undef Trans
#undef field_Field
} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_SpectralFunctionSpace_h
