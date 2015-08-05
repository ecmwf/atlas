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
 * \brief Nodes class that owns a collection of fields
 */
/// TEMPORARILY DERIVE FROM FunctionSpace UNTIL NEXT DESIGN IS COMPLETE
class Nodes : public FunctionSpace {

public: // methods

//-- Constructors

  Nodes(size_t size);

//-- Accessors

  virtual const Field& field(const std::string& name) const;
  virtual       Field& field(const std::string& name);
  virtual bool has_field(const std::string& name) const { return (fields_.find(name) != fields_.end()); }

  virtual const Field& field(size_t) const;
  virtual       Field& field(size_t);
  virtual size_t nb_fields() const { return fields_.size(); }

//  const Metadata& metadata() const { return metadata_; }
//        Metadata& metadata()       { return metadata_; }

  const Field& global_index() const { return *global_index_; }
        Field& global_index()       { return *global_index_; }

  const Field& remote_index() const { return *remote_index_; }
        Field& remote_index()       { return *remote_index_; }

  const Field& partition() const { return *partition_; }
        Field& partition()       { return *partition_; }

  const Field& ghost() const { return *ghost_; }
        Field& ghost()       { return *ghost_; }

  const Field& halo() const { return *halo_; }
        Field& halo()       { return *halo_; }

  const Field& topology() const { return *topology_; }
        Field& topology()       { return *topology_; }

  const Field& lonlat() const { return *lonlat_; }
        Field& lonlat()       { return *lonlat_; }

  size_t size() const { return dof_; }

// -- Modifiers

  virtual Field& add( Field* ); // Take ownership!

  void resize( size_t );

private:

  void remove_field(const std::string& name);

  template< typename DATA_TYPE >
  Field& create_field(const std::string& name, size_t nb_vars, CreateBehavior b = IF_EXISTS_FAIL );

  const std::string& name() const;

  int index() const;

  // This is a Fortran view of the shape (i.e. reverse order)
  const std::vector<int>& shapef() const;

  const std::vector<size_t>& shape() const;
  size_t shape(const size_t i) const;
  void resize( const std::vector<size_t>& shape );

  const Mesh& mesh() const;
  Mesh& mesh();

  size_t dof() const;

  size_t glb_dof() const;

  void print(std::ostream&, bool dump = false) const;

private:

  typedef std::map< std::string, eckit::SharedPtr<Field> >  FieldMap;

private:

  //size_t size_;
  FieldMap fields_;
  //Metadata metadata_;

  // Cached shortcuts to specific fields in fields_
  Field* global_index_;
  Field* remote_index_;
  Field* partition_;
  Field* ghost_;
  Field* halo_;
  Field* topology_;
  Field* lonlat_;

};

//------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif
