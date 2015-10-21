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
/// @date August 2015

#ifndef atlas_Nodes_H
#define atlas_Nodes_H

#include <string>
#include <map>
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/Metadata.h"
#include "atlas/FunctionSpace.h"

namespace atlas { class Field; }

namespace atlas {

/**
 * \brief Nodes class that owns a collection of fields defined in nodes of the mesh
 */
class Nodes : public eckit::Owned {

public: // methods

//-- Constructors

  /// @brief Construct "size" nodes
  Nodes(size_t size);

//-- Accessors

  virtual const Field& field(const std::string& name) const;
  virtual       Field& field(const std::string& name);
  virtual bool has_field(const std::string& name) const { return (fields_.find(name) != fields_.end()); }

  virtual const Field& field(size_t) const;
  virtual       Field& field(size_t);
  virtual size_t nb_fields() const { return fields_.size(); }

  const Metadata& metadata() const { return metadata_; }
        Metadata& metadata()       { return metadata_; }

  const Field& global_index() const { return *global_index_; }
        Field& global_index()       { return *global_index_; }

  const Field& remote_index() const { return *remote_index_; }
        Field& remote_index()       { return *remote_index_; }

  const Field& partition() const { return *partition_; }
        Field& partition()       { return *partition_; }

  const Field& lonlat() const { return *lonlat_; }
        Field& lonlat()       { return *lonlat_; }

  const Field& ghost() const { return *ghost_; }
        Field& ghost()       { return *ghost_; }

  size_t size() const { return size_; }

// -- Modifiers

  virtual Field& add( Field* ); // Take ownership!

  void resize( size_t );

  void remove_field(const std::string& name);

private:

  void print(std::ostream&) const;

  friend std::ostream& operator<<(std::ostream& s, const Nodes& p) {
    p.print(s);
    return s;
  }

private:

  typedef std::map< std::string, eckit::SharedPtr<Field> >  FieldMap;

private:

  size_t size_;
  FieldMap fields_;
  Metadata metadata_;

  // Cached shortcuts to specific fields in fields_
  Field* global_index_;
  Field* remote_index_;
  Field* partition_;
  Field* lonlat_;
  Field* ghost_;

};

#define Char char
extern "C"
{
int atlas__Nodes__size (Nodes* This);
void atlas__Nodes__resize (Nodes* This, int size);
int atlas__Nodes__nb_fields (Nodes* This);
void atlas__Nodes__add (Nodes* This, Field* field);
void atlas__Nodes__remove_field (Nodes* This, char* name);
int atlas__Nodes__has_field (Nodes* This, char* name);
Field* atlas__Nodes__field_by_name (Nodes* This, char* name);
Field* atlas__Nodes__field_by_idx (Nodes* This, int idx);
Metadata* atlas__Nodes__metadata(Nodes* This);
void atlas__Nodes__str (Nodes* This, Char* &str, int &size);
}
#undef Char

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
