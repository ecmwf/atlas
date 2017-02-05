/*
 * (C) Copyright 1996-2017 ECMWF.
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

#include <map>
#include <string>
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/util/Metadata.h"

namespace atlas { namespace field { class Field; } }

namespace atlas {
namespace mesh {

/**
 * \brief Nodes class that owns a collection of fields defined in nodes of the mesh
 */
class Nodes : public eckit::Owned {
public:
  typedef IrregularConnectivity Connectivity;

public: // methods

//-- Constructors

  /// @brief Construct "size" nodes
  Nodes();
//  Nodes(size_t size);

//-- Accessors

  const field::Field& field(const std::string& name) const;
        field::Field& field(const std::string& name);
  bool has_field(const std::string& name) const { return (fields_.find(name) != fields_.end()); }

  const field::Field& field(size_t) const;
        field::Field& field(size_t);
  size_t nb_fields() const { return fields_.size(); }

  const util::Metadata& metadata() const { return metadata_; }
        util::Metadata& metadata()       { return metadata_; }

  const field::Field& global_index() const { return *global_index_; }
        field::Field& global_index()       { return *global_index_; }

  const field::Field& remote_index() const { return *remote_index_; }
        field::Field& remote_index()       { return *remote_index_; }

  const field::Field& partition() const { return *partition_; }
        field::Field& partition()       { return *partition_; }

  const field::Field& lonlat() const { return *lonlat_; }
        field::Field& lonlat()       { return *lonlat_; }

  const field::Field& ghost() const { return *ghost_; }
        field::Field& ghost()       { return *ghost_; }

  /// @brief Node to Edge connectivity table
  const Connectivity& edge_connectivity() const;
        Connectivity& edge_connectivity();

  /// @brief Node to Cell connectivity table
  const Connectivity& cell_connectivity() const;
        Connectivity& cell_connectivity();

  const Connectivity& connectivity(const std::string& name) const;
        Connectivity& connectivity(const std::string& name);

  size_t size() const { return size_; }

// -- Modifiers

  field::Field& add( field::Field* ); // Take ownership!

  void resize( size_t );

  void remove_field(const std::string& name);

  Connectivity& add( mesh::Connectivity* );

  /// @brief Return the memory footprint of the Nodes
  size_t footprint() const;

private:

  void print(std::ostream&) const;

  friend std::ostream& operator<<(std::ostream& s, const Nodes& p) {
    p.print(s);
    return s;
  }

private:

  typedef std::map< std::string, eckit::SharedPtr<field::Field> >  FieldMap;
  typedef std::map< std::string, eckit::SharedPtr<Connectivity> >  ConnectivityMap;

private:

  size_t size_;
  FieldMap fields_;
  ConnectivityMap connectivities_;

  util::Metadata metadata_;

  // Cached shortcuts to specific fields in fields_
  field::Field* global_index_;
  field::Field* remote_index_;
  field::Field* partition_;
  field::Field* lonlat_;
  field::Field* ghost_;

// Cached shortcuts to specific connectivities in connectivities_
  Connectivity* edge_connectivity_;
  Connectivity* cell_connectivity_;

};

inline const Nodes::Connectivity& Nodes::edge_connectivity() const
{
  return *edge_connectivity_;
}

inline Nodes::Connectivity& Nodes::edge_connectivity()
{
  return *edge_connectivity_;
}

inline const Nodes::Connectivity& Nodes::cell_connectivity() const
{
  return *cell_connectivity_;
}

inline Nodes::Connectivity& Nodes::cell_connectivity()
{
  return *cell_connectivity_;
}

extern "C"
{
Nodes* atlas__mesh__Nodes__create();
void atlas__mesh__Nodes__delete (Nodes* This);
size_t atlas__mesh__Nodes__size (Nodes* This);
void atlas__mesh__Nodes__resize (Nodes* This, size_t size);
size_t atlas__mesh__Nodes__nb_fields (Nodes* This);
void atlas__mesh__Nodes__add_field (Nodes* This, field::Field* field);
void atlas__mesh__Nodes__remove_field (Nodes* This, char* name);
int atlas__mesh__Nodes__has_field (Nodes* This, char* name);
field::Field* atlas__mesh__Nodes__field_by_name (Nodes* This, char* name);
field::Field* atlas__mesh__Nodes__field_by_idx (Nodes* This, size_t idx);
util::Metadata* atlas__mesh__Nodes__metadata(Nodes* This);
void atlas__mesh__Nodes__str (Nodes* This, char* &str, int &size);
IrregularConnectivity* atlas__mesh__Nodes__edge_connectivity(Nodes* This);
IrregularConnectivity* atlas__mesh__Nodes__cell_connectivity(Nodes* This);
IrregularConnectivity* atlas__mesh__Nodes__connectivity (Nodes* This, char* name);
void atlas__mesh__Nodes__add_connectivity (Nodes* This, IrregularConnectivity* connectivity);
field::Field* atlas__mesh__Nodes__lonlat(Nodes* This);
field::Field* atlas__mesh__Nodes__global_index(Nodes* This);
field::Field* atlas__mesh__Nodes__remote_index(Nodes* This);
field::Field* atlas__mesh__Nodes__partition(Nodes* This);
field::Field* atlas__mesh__Nodes__ghost(Nodes* This);
}

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif
