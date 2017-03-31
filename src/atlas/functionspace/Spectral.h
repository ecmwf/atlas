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

#include "atlas/library/config.h"
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
namespace detail {

// -------------------------------------------------------------------

class Spectral : public FunctionSpaceImpl
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

  void norm( const field::Field&, double& norm, int rank=0 ) const;
  void norm( const field::Field&, double norm_per_level[], int rank=0 ) const;
  void norm( const field::Field&, std::vector<double>& norm_per_level, int rank=0 ) const;

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

} // detail

// -------------------------------------------------------------------

class Spectral : public FunctionSpace {
public:

  Spectral( const FunctionSpace& );
  Spectral( const size_t truncation );
  Spectral( trans::Trans& );


  operator bool() const { return valid(); }
  bool valid() const { return functionspace_; }

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

  void norm( const field::Field&, double& norm, int rank=0 ) const;
  void norm( const field::Field&, double norm_per_level[], int rank=0 ) const;
  void norm( const field::Field&, std::vector<double>& norm_per_level, int rank=0 ) const;

  size_t nb_spectral_coefficients() const;
  size_t nb_spectral_coefficients_global() const;

private:

  const detail::Spectral* functionspace_;
};


// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  const detail::Spectral* atlas__SpectralFunctionSpace__new__truncation (int truncation);
  const detail::Spectral* atlas__SpectralFunctionSpace__new__trans (trans::Trans* trans);
  void atlas__SpectralFunctionSpace__delete (detail::Spectral* This);
  field::Field* atlas__fs__Spectral__create_field_name_kind(const detail::Spectral* This, const char* name, int kind, const eckit::Parametrisation* options);
  field::Field* atlas__fs__Spectral__create_field_name_kind_lev(const detail::Spectral* This, const char* name, int kind, int levels, const eckit::Parametrisation* options);
  void atlas__SpectralFunctionSpace__gather(const detail::Spectral* This, const field::Field* local, field::Field* global);
  void atlas__SpectralFunctionSpace__gather_fieldset(const detail::Spectral* This, const field::FieldSet* local, field::FieldSet* global);
  void atlas__SpectralFunctionSpace__scatter(const detail::Spectral* This, const field::Field* global, field::Field* local);
  void atlas__SpectralFunctionSpace__scatter_fieldset(const detail::Spectral* This, const field::FieldSet* global, field::FieldSet* local);
  void atlas__SpectralFunctionSpace__norm(const detail::Spectral* This, const field::Field* field, double norm[], int rank);
}

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_SpectralFunctionSpace_h
