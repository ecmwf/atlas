/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_SpectralFunctionSpace_h
#define atlas_functionspace_SpectralFunctionSpace_h

#include "atlas/internals/atlas_defines.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace field {
    class Field;
    class FieldSet;
}
}

namespace atlas {
namespace trans {
    class Trans;
}
}

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class Spectral : public FunctionSpace
{
public:

    typedef eckit::SharedPtr<Spectral> Ptr;

public:

  Spectral(const size_t truncation);

  Spectral(trans::Trans& );

  virtual ~Spectral();

  virtual std::string name() const { return "Spectral"; }

  /// @brief Create a spectral field
  template <typename DATATYPE> field::Field* createField(
    const std::string& name,
    const eckit::Parametrisation& = util::NoConfig() ) const;
  template <typename DATATYPE> field::Field* createField(
    const std::string& name,
    size_t levels,
    const eckit::Parametrisation& = util::NoConfig() ) const;

  void gather( const field::FieldSet&, field::FieldSet& ) const;
  void gather( const field::Field&,    field::Field& ) const;

  void scatter( const field::FieldSet&, field::FieldSet& ) const;
  void scatter( const field::Field&,    field::Field& ) const;

  std::string checksum( const field::FieldSet& ) const;
  std::string checksum( const field::Field& ) const;

public: // methods

  size_t nb_spectral_coefficients() const;
  size_t nb_spectral_coefficients_global() const;

private: // methods

  size_t config_size(const eckit::Parametrisation& config) const;
  size_t footprint() const;

private: // data

  size_t truncation_;

  trans::Trans* trans_;

};

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
#define Trans trans::Trans
#define field_Field field::Field
#define field_FieldSet field::FieldSet
#define Options eckit::Parametrisation
extern "C"
{
  Spectral* atlas__SpectralFunctionSpace__new__truncation (int truncation);
  Spectral* atlas__SpectralFunctionSpace__new__trans (Trans* trans);
  void atlas__SpectralFunctionSpace__delete (Spectral* This);
  field_Field* atlas__fs__Spectral__create_field_name_kind(const Spectral* This, const char* name, int kind, const Options* options);
  field_Field* atlas__fs__Spectral__create_field_name_kind_lev(const Spectral* This, const char* name, int kind, int levels, const Options* options);
  void atlas__SpectralFunctionSpace__gather(const Spectral* This, const field_Field* local, field_Field* global);
  void atlas__SpectralFunctionSpace__gather_fieldset(const Spectral* This, const field_FieldSet* local, field_FieldSet* global);
  void atlas__SpectralFunctionSpace__scatter(const Spectral* This, const field_Field* global, field_Field* local);
  void atlas__SpectralFunctionSpace__scatter_fieldset(const Spectral* This, const field_FieldSet* global, field_FieldSet* local);
}
#undef Trans
#undef field_Field
#undef field_FieldSet
#undef Options
} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_SpectralFunctionSpace_h
