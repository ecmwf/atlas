/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_functionspace__StructuredColumns_h
#define atlas_functionspace_functionspace__StructuredColumns_h

#include "atlas/internals/atlas_defines.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/util/Config.h"
#include "atlas/grid/Grid.h"
#include "atlas/array/DataType.h"

namespace atlas {
namespace parallel {
    class GatherScatter;
    class Checksum;
}
}

namespace atlas {
namespace field {
    class Field;
    class FieldSet;
}
}

namespace atlas {
namespace grid {
    class Structured;
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

class StructuredColumns : public FunctionSpace
{
public:

    typedef eckit::SharedPtr<StructuredColumns> Ptr;

public:

  StructuredColumns( const grid::Grid& );

  virtual ~StructuredColumns();

  virtual std::string name() const { return "StructuredColumns"; }

  /// @brief Create a Structured field
  field::Field* createField(
      const std::string& name,
      array::DataType,
      const eckit::Parametrisation& = util::NoConfig() ) const;
  field::Field* createField(
      const std::string& name,
      array::DataType,
      size_t levels,
      const eckit::Parametrisation& = util::NoConfig() ) const;
  template <typename DATATYPE> field::Field* createField(
      const std::string& name,
      const eckit::Parametrisation& = util::NoConfig() ) const;
  template <typename DATATYPE> field::Field* createField(
      const std::string& name,
      size_t levels,
      const eckit::Parametrisation& = util::NoConfig() ) const;

  void gather( const field::FieldSet&, field::FieldSet& ) const;
  void gather( const field::Field&, field::Field& ) const;

  void scatter( const field::FieldSet&, field::FieldSet& ) const;
  void scatter( const field::Field&, field::Field& ) const;

  size_t npts() const            { return npts_; }
  size_t nlat() const            { return nlat_; }
  size_t nlon(size_t jlat) const { return nlon_[jlat]; }

  double lat(size_t jlat) const;
  double lon(size_t jlat, size_t jlon) const;

  std::string checksum( const field::FieldSet& ) const;
  std::string checksum( const field::Field& ) const;

  const grid::Structured& grid() const { return grid_; }

private: // methods

  size_t config_size(const eckit::Parametrisation& config) const;

private: // data

  size_t npts_;
  size_t nlat_;
  size_t first_lat_;
  std::vector<size_t> nlon_;
  std::vector<size_t> first_lon_;

  trans::Trans* trans_;
  const grid::Structured grid_;
  parallel::GatherScatter* gather_scatter_;
  parallel::Checksum* checksum_;

};

// -------------------------------------------------------------------
// inline methods

template <typename DATATYPE>
inline field::Field* StructuredColumns::createField(
    const std::string& name,
    const eckit::Parametrisation& options) const
{
  return createField(name,array::DataType::create<DATATYPE>(),options);
}

template <typename DATATYPE>
inline field::Field* StructuredColumns::createField(
    const std::string& name,
    size_t levels,
    const eckit::Parametrisation& options) const
{
  return createField(name,array::DataType::create<DATATYPE>(),levels,options);
}

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
#define Char char
#define grid_Grid grid::Grid::grid_t
#define field_Field field::Field
#define field_FieldSet field::FieldSet
#define Options eckit::Parametrisation
extern "C"
{
  StructuredColumns* atlas__functionspace__StructuredColumns__new__grid (const grid_Grid* grid);
  void atlas__functionspace__StructuredColumns__delete (StructuredColumns* This);
  field_Field* atlas__fs__StructuredColumns__create_field_name_kind (const StructuredColumns* This, const char* name, int kind, const Options* options);
  field_Field* atlas__fs__StructuredColumns__create_field_name_kind_lev (const StructuredColumns* This, const char* name, int kind, int levels, const Options* options);
  void atlas__functionspace__StructuredColumns__gather (const StructuredColumns* This, const field_Field* local, field_Field* global);
  void atlas__functionspace__StructuredColumns__scatter (const StructuredColumns* This, const field_Field* global, field_Field* local);
  void atlas__fs__StructuredColumns__checksum_fieldset(const StructuredColumns* This, const field_FieldSet* fieldset, Char* &checksum, int &size, int &allocated);
  void atlas__fs__StructuredColumns__checksum_field(const StructuredColumns* This, const field_Field* field, Char* &checksum, int &size, int &allocated);
}

#undef grid_Grid
#undef field_FieldSet
#undef field_Field
#undef Char
#undef Options


} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_functionspace__StructuredColumns_h
