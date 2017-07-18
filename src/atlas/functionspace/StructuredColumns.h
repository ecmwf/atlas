/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_functionspace__StructuredColumns_h
#define atlas_functionspace_functionspace__StructuredColumns_h

#include "atlas/library/config.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/util/Config.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/array/DataType.h"
#include "atlas/field/Field.h"

namespace atlas {
namespace parallel {
    class GatherScatter;
    class Checksum;
}
}

namespace atlas {
  class FieldSet;
namespace field {
    class FieldSetImpl;
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

namespace detail {
class StructuredColumns : public FunctionSpaceImpl {

public:

  StructuredColumns( const Grid& );

  StructuredColumns( const Grid&, const grid::Partitioner& );

  virtual ~StructuredColumns();

  virtual std::string name() const { return "StructuredColumns"; }

  /// @brief Create a Structured field
  Field createField(
      const std::string& name,
      array::DataType,
      const eckit::Parametrisation& = util::NoConfig() ) const;
  Field createField(
      const std::string& name,
      array::DataType,
      size_t levels,
      const eckit::Parametrisation& = util::NoConfig() ) const;
  template <typename DATATYPE> Field createField(
      const std::string& name,
      const eckit::Parametrisation& = util::NoConfig() ) const;
  template <typename DATATYPE> Field createField(
      const std::string& name,
      size_t levels,
      const eckit::Parametrisation& = util::NoConfig() ) const;

  void gather( const FieldSet&, FieldSet& ) const;
  void gather( const Field&, Field& ) const;

  void scatter( const FieldSet&, FieldSet& ) const;
  void scatter( const Field&, Field& ) const;

  size_t size() const       { return npts_; }
  size_t ny() const         { return ny_; }
  size_t nx(size_t j) const { return nx_[j]; }

  double x(size_t i, size_t j) const;
  double y(size_t j) const;

  std::string checksum( const FieldSet& ) const;
  std::string checksum( const Field& ) const;

  const grid::StructuredGrid& grid() const { return grid_; }

private: // methods

  size_t config_size(const eckit::Parametrisation& config) const;
  size_t footprint() const;

private: // data

  size_t npts_;
  size_t ny_;
  size_t j_begin_;
  std::vector<size_t> nx_;
  std::vector<size_t> i_begin_;

  trans::Trans* trans_;
  const grid::StructuredGrid grid_;
  parallel::GatherScatter* gather_scatter_;
  parallel::Checksum* checksum_;

};

// -------------------------------------------------------------------
// inline methods

template <typename DATATYPE>
inline Field StructuredColumns::createField(
    const std::string& name,
    const eckit::Parametrisation& options) const
{
  return createField(name,array::DataType::create<DATATYPE>(),options);
}

template <typename DATATYPE>
inline Field StructuredColumns::createField(
    const std::string& name,
    size_t levels,
    const eckit::Parametrisation& options) const
{
  return createField(name,array::DataType::create<DATATYPE>(),levels,options);
}
}

// -------------------------------------------------------------------

class StructuredColumns: public FunctionSpace {

public:

  StructuredColumns();
  StructuredColumns( const FunctionSpace& );
  StructuredColumns( const Grid& );
  StructuredColumns( const Grid&, const grid::Partitioner& );

  operator bool() const { return valid(); }
  bool valid() const { return functionspace_; }

  size_t size() const            { return functionspace_->size(); }
  size_t ny() const            { return functionspace_->ny(); }
  size_t nx(size_t j) const { return functionspace_->nx(j); }

  double y(size_t j) const { return functionspace_->y(j); }
  double x(size_t i, size_t j) const { return functionspace_->x(i,j); }

  const grid::StructuredGrid& grid() const { return functionspace_->grid(); }

  Field createField(
      const std::string& name,
      array::DataType,
      const eckit::Parametrisation& = util::NoConfig() ) const;
  Field createField(
      const std::string& name,
      array::DataType,
      size_t levels,
      const eckit::Parametrisation& = util::NoConfig() ) const;
  template <typename DATATYPE> Field createField(
      const std::string& name,
      const eckit::Parametrisation& = util::NoConfig() ) const;
  template <typename DATATYPE> Field createField(
      const std::string& name,
      size_t levels,
      const eckit::Parametrisation& = util::NoConfig() ) const;

  void gather( const FieldSet&, FieldSet& ) const;
  void gather( const Field&, Field& ) const;

  void scatter( const FieldSet&, FieldSet& ) const;
  void scatter( const Field&, Field& ) const;

  std::string checksum( const FieldSet& ) const;
  std::string checksum( const Field& ) const;

private:

  const detail::StructuredColumns* functionspace_;
};

// -------------------------------------------------------------------
// inline methods

template <typename DATATYPE>
inline Field StructuredColumns::createField(
    const std::string& name,
    const eckit::Parametrisation& options) const
{
  return functionspace_->createField<DATATYPE>(name,options);
}

template <typename DATATYPE>
inline Field StructuredColumns::createField(
    const std::string& name,
    size_t levels,
    const eckit::Parametrisation& options) const
{
  return functionspace_->createField<DATATYPE>(name,levels,options);
}

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C"
{
  const detail::StructuredColumns* atlas__functionspace__StructuredColumns__new__grid (const Grid::Implementation* grid);
  void atlas__functionspace__StructuredColumns__delete (detail::StructuredColumns* This);
  field::FieldImpl* atlas__fs__StructuredColumns__create_field_name_kind (const detail::StructuredColumns* This, const char* name, int kind, const eckit::Parametrisation* options);
  field::FieldImpl* atlas__fs__StructuredColumns__create_field_name_kind_lev (const detail::StructuredColumns* This, const char* name, int kind, int levels, const eckit::Parametrisation* options);
  void atlas__functionspace__StructuredColumns__gather (const detail::StructuredColumns* This, const field::FieldImpl* local, field::FieldImpl* global);
  void atlas__functionspace__StructuredColumns__scatter (const detail::StructuredColumns* This, const field::FieldImpl* global, field::FieldImpl* local);
  void atlas__fs__StructuredColumns__checksum_fieldset(const detail::StructuredColumns* This, const field::FieldSetImpl* fieldset, char* &checksum, int &size, int &allocated);
  void atlas__fs__StructuredColumns__checksum_field(const detail::StructuredColumns* This, const field::FieldImpl* field, char* &checksum, int &size, int &allocated);
}

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_functionspace__StructuredColumns_h
