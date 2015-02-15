/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @author Peter Bispham
/// @author Tiago Quintino
/// @author Pedro Maciel
/// @date Jan 2015

#ifndef atlas_FieldSet_H
#define atlas_FieldSet_H

#include <vector>

#include "atlas/atlas_config.h"

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "eckit/types/Types.h"

#include "atlas/Field.h"
#include "atlas/GridSet.h"
#include "atlas/Grid.h"
#include "atlas/Mesh.h"
#include "atlas/Metadata.h"


// forward declarations
namespace eckit {
  class PathName;
  class DataHandle;
  namespace grib { class GribHandle; }
}

namespace atlas {

class FieldSet;
class Field;
std::vector<Field*>& __private_get_raw_fields_ptr (FieldSet* This);

/**
 * @brief Represents a set of fields, where order is preserved (no ownership)
 */
class FieldSet : public eckit::Owned {

public: // types

  typedef eckit::SharedPtr< FieldSet > Ptr;

public: // methods

  /// Constructs an empty FieldSet
  FieldSet(const std::string& name = "untitled");

  /// Constructs from a path (e.g. a GRIB file)
  FieldSet(const eckit::PathName&);

  /// Constructs from a buffer
  FieldSet(const eckit::Buffer&);

  /// Constructs from a DataHandle
  FieldSet(const eckit::DataHandle&);

  /// Constructs from predefined fields (takes ownership of the fields)
  FieldSet(const Field::Vector& fields);

  size_t size() const { return  fields_.size(); }
  bool empty()  const { return !fields_.size(); }

  const std::string& name() const { return name_; }
        std::string& name()       { return name_; }

  const Field& operator[](const size_t& i) const { return field(i); }
        Field& operator[](const size_t& i)       { return field(i); }

  const Field& field(const size_t& i) const { ASSERT(i<size()); return *fields_[i]; }
        Field& field(const size_t& i)       { ASSERT(i<size()); return *fields_[i]; }

  const Field::Vector& fields() const { return fields_; }
        Field::Vector& fields()       { return fields_; }

  std::vector< std::string > field_names() const;

  void add_field(Field::Ptr field);

  bool has_field(const std::string& name) const;

  Field& field(const std::string& name) const;
  
private: // methods

  /// @todo to be removed
  Field::Ptr create_field(eckit::grib::GribHandle&);

  /// Check if all fields share same grid
  bool haveSameGrid() const;

protected: // data

  atlas::GridSet                  gridset_; ///< @todo remove this gridset ownership into the mir::FieldSet

  Field::Vector                   fields_;  ///< field handle storage
  std::string                     name_;    ///< internal name
  std::map< std::string, size_t > index_;   ///< name-to-index map, to refer fields by name

private:

  // In order to return raw pointers to C interface
  friend std::vector<Field*>& __private_get_raw_fields_ptr (FieldSet* This);
  std::vector<Field*> fields_raw_ptr_;

};


// C wrapper interfaces to C++ routines
extern "C"
{
  FieldSet* atlas__FieldSet__new           (char* name);
  void      atlas__FieldSet__delete        (FieldSet* This);
  void      atlas__FieldSet__fields        (FieldSet* This, Field** &fields, int &nb_fields);
  void      atlas__FieldSet__add_field     (FieldSet* This, Field* field);
  int       atlas__FieldSet__has_field     (FieldSet* This, char* name);
  int       atlas__FieldSet__size          (FieldSet* This);
  Field*    atlas__FieldSet__field_by_name (FieldSet* This, char* name);
  Field*    atlas__FieldSet__field_by_idx  (FieldSet* This, int idx);
}


} // namespace atlas


#endif
