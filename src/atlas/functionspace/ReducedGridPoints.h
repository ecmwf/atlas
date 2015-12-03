/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_functionspace_functionspace__ReducedGridPoints_h
#define atlas_functionspace_functionspace__ReducedGridPoints_h

#include "atlas/atlas_defines.h"
#include "atlas/FunctionSpace.h"

namespace atlas { namespace trans { class Trans; } }
namespace atlas { namespace mpl { class GatherScatter; } }
namespace atlas { namespace mpl { class Checksum; } }
namespace atlas { class Field; }
namespace atlas { class FieldSet; }

namespace atlas {
namespace functionspace {

// -------------------------------------------------------------------

class ReducedGridPoints : public next::FunctionSpace
{
public:

  ReducedGridPoints( const Grid& );

  virtual ~ReducedGridPoints();

  virtual std::string name() const { return "ReducedGrid"; }

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

  std::string checksum( const FieldSet& ) const;
  std::string checksum( const Field& ) const;

private: // data

  trans::Trans* trans_;
  const Grid& grid_;
  mpl::GatherScatter* gather_scatter_;
  mpl::Checksum* checksum_;

};

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
#define Char char
extern "C"
{
  ReducedGridPoints* atlas__functionspace__ReducedGridPoints__new__grid (const Grid* grid);
  void atlas__functionspace__ReducedGridPoints__delete (ReducedGridPoints* This);

  Field* atlas__functionspace__ReducedGridPoints__create_field (const ReducedGridPoints* This, const char* name);

  Field* atlas__functionspace__ReducedGridPoints__create_field_lev (const ReducedGridPoints* This, const char* name, int levels);

  Field* atlas__functionspace__ReducedGridPoints__create_gfield (const ReducedGridPoints* This, const char* name);

  Field* atlas__functionspace__ReducedGridPoints__create_gfield_lev (const ReducedGridPoints* This, const char* name, int levels);

  void atlas__functionspace__ReducedGridPoints__gather (const ReducedGridPoints* This, const Field* local, Field* global);

  void atlas__functionspace__ReducedGridPoints__scatter (const ReducedGridPoints* This, const Field* global, Field* local);
  
  void atlas__fs__ReducedGridPoints__checksum_fieldset(const ReducedGridPoints* This, const FieldSet* fieldset, Char* &checksum, int &size, int &allocated);

  void atlas__fs__ReducedGridPoints__checksum_field(const ReducedGridPoints* This, const Field* field, Char* &checksum, int &size, int &allocated);
  

}
#undef Char

} // namespace functionspace
} // namespace atlas

#endif // atlas_functionspace_functionspace__ReducedGridPoints_h
