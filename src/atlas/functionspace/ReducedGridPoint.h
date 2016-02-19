/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_functionspace__ReducedGridPoint_h
#define atlas_functionspace_functionspace__ReducedGridPoint_h

#include "atlas/atlas_defines.h"
#include "atlas/functionspace/FunctionSpace.h"

namespace atlas { namespace trans { class Trans; } }
namespace atlas { namespace mpl { class GatherScatter; } }
namespace atlas { namespace mpl { class Checksum; } }
namespace atlas { namespace grids { class ReducedGrid; } }
namespace atlas { class Field; }
namespace atlas { class FieldSet; }

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class ReducedGridPoint : public FunctionSpace
{
public:

  ReducedGridPoint( const Grid& );

  virtual ~ReducedGridPoint();

  virtual std::string name() const { return "ReducedGridPoint"; }

  /// @brief Create a ReducedGrid field
  template <typename DATATYPE> Field* createField(const std::string& name) const;
  template <typename DATATYPE> Field* createField(const std::string& name, size_t levels) const;

  /// @brief Create a global ReducedGrid field
  template <typename DATATYPE> Field* createGlobalField(const std::string& name) const;
  template <typename DATATYPE> Field* createGlobalField(const std::string& name, size_t levels) const;

  void gather( const FieldSet&, FieldSet& ) const;
  void gather( const Field&, Field& ) const;

  void scatter( const FieldSet&, FieldSet& ) const;
  void scatter( const Field&, Field& ) const;

  size_t npts() const            { return npts_; }
  size_t nlat() const            { return nlat_; }
  size_t nlon(size_t jlat) const { return nlon_[jlat]; }

  double lat(size_t jlat) const;
  double lon(size_t jlat, size_t jlon) const;

  std::string checksum( const FieldSet& ) const;
  std::string checksum( const Field& ) const;

  const grids::ReducedGrid& grid() const { return *grid_; }

private: // data

  size_t npts_;
  size_t nlat_;
  size_t first_lat_;
  std::vector<size_t> nlon_;
  std::vector<size_t> first_lon_;

  trans::Trans* trans_;
  const grids::ReducedGrid* grid_;
  mpl::GatherScatter* gather_scatter_;
  mpl::Checksum* checksum_;

};

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
#define Char char
extern "C"
{
  ReducedGridPoint* atlas__functionspace__ReducedGridPoint__new__grid (const Grid* grid);
  void atlas__functionspace__ReducedGridPoint__delete (ReducedGridPoint* This);

  Field* atlas__functionspace__ReducedGridPoint__create_field (const ReducedGridPoint* This, const char* name);

  Field* atlas__functionspace__ReducedGridPoint__create_field_lev (const ReducedGridPoint* This, const char* name, int levels);

  Field* atlas__functionspace__ReducedGridPoint__create_gfield (const ReducedGridPoint* This, const char* name);

  Field* atlas__functionspace__ReducedGridPoint__create_gfield_lev (const ReducedGridPoint* This, const char* name, int levels);

  void atlas__functionspace__ReducedGridPoint__gather (const ReducedGridPoint* This, const Field* local, Field* global);

  void atlas__functionspace__ReducedGridPoint__scatter (const ReducedGridPoint* This, const Field* global, Field* local);
  
  void atlas__fs__ReducedGridPoint__checksum_fieldset(const ReducedGridPoint* This, const FieldSet* fieldset, Char* &checksum, int &size, int &allocated);

  void atlas__fs__ReducedGridPoint__checksum_field(const ReducedGridPoint* This, const Field* field, Char* &checksum, int &size, int &allocated);
  

}
#undef Char

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_functionspace__ReducedGridPoint_h
