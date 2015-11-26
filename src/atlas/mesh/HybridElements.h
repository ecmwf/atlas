/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @file HybridElements.h
/// @author Willem Deconinck
/// @date October 2015
///
/// This file describes the HybridElements class for a Mesh.

#ifndef atlas_mesh_HybridElements_H
#define atlas_mesh_HybridElements_H

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/Connectivity.h"
#include "atlas/Metadata.h"

namespace atlas { namespace mesh { class ElementType; } }
namespace atlas { namespace mesh { class Elements; } }
namespace atlas { class Field; }

namespace atlas {
namespace mesh {

// -------------------------------------------------------------------------------

/// @brief HybridElements class that describes elements of different types
class HybridElements : public eckit::Owned {
friend class Elements;
public:
  typedef MultiBlockConnectivity Connectivity;

public: // methods

//-- Constructors

  HybridElements();
  virtual ~HybridElements();

//-- Accessors

  /// @brief Number of elements
  size_t size() const;

  /// @brief Number of nodes for given element
  size_t nb_nodes( size_t elem_idx ) const;

  /// @brief Number of edges for given element
  size_t nb_edges( size_t elem_idx ) const;

  /// @brief Element type index for given element
  size_t type_idx( size_t elem_idx ) const;

  /// @brief Element type name for given element
  const std::string& name( size_t elem_idx ) const;

  /// @brief Element to Node connectivity table
  const HybridElements::Connectivity& node_connectivity() const;
        HybridElements::Connectivity& node_connectivity();

  /// @brief Element to Edge connectivity table
  const HybridElements::Connectivity& edge_connectivity() const;
        HybridElements::Connectivity& edge_connectivity();

  /// @brief Element to Cell connectivity table
  const HybridElements::Connectivity& cell_connectivity() const;
        HybridElements::Connectivity& cell_connectivity();

  /// @brief Number of types present in HybridElements
  size_t nb_types() const;

  /// @brief The element_type description for given type
  const ElementType& element_type( size_t type_idx ) const;

  /// @brief Sub-elements convenience class for given type
  /// This allows optimized access to connectivities and loops.
  const Elements& elements( size_t type_idx ) const;
        Elements& elements( size_t type_idx );

// -- Modifiers

  /// @brief Add a new element type with given number of elements and node-connectivity
  /// @return type_idx of the added element type
  size_t add( const ElementType*, size_t nb_elements, const std::vector<idx_t> &node_connectivity );

  /// @brief Add a new element type with given number of elements and node-connectivity
  /// @return type_idx of the added element type
  size_t add( const ElementType*, size_t nb_elements, const idx_t node_connectivity[] );

  /// @brief Add a new element type with given number of elements and node-connectivity
  /// @return type_idx of the added element type
  size_t add( const ElementType*, size_t nb_elements, const idx_t node_connectivity[], bool fortran_array );

  /// @brief Add a new element type from existing Elements.
  /// Data will be copied.
  /// @return type_idx of the added element type
  size_t add( const Elements& );

private:

// -- Data
  size_t size_;                          //!< total number of elements

// -- Data: one value per type
  std::vector<size_t> elements_size_;
  std::vector<size_t> elements_begin_;
  std::vector< eckit::SharedPtr<const ElementType> > element_types_;

// -- Data: one value per element
  std::vector<size_t> type_idx_;

// -- Connectivity tables
  Connectivity* node_connectivity_;
  Connectivity* edge_connectivity_;
  Connectivity* cell_connectivity_;
  
// -- Sub elements
  std::vector< eckit::SharedPtr<Elements> > elements_;

// -- New stuff

private:

  typedef std::map< std::string, eckit::SharedPtr<Field>        >  FieldMap;
  typedef std::map< std::string, eckit::SharedPtr<Connectivity> >  ConnectivityMap;

  void resize( size_t size );

public:
  Field& add( Field* field );
  void remove_field(const std::string& name);

  const Field& field(const std::string& name) const;
        Field& field(const std::string& name);
  bool has_field(const std::string& name) const { return (fields_.find(name) != fields_.end()); }

  const Field& field(size_t) const;
        Field& field(size_t);
  size_t nb_fields() const { return fields_.size(); }

  const Metadata& metadata() const { return metadata_; }
        Metadata& metadata()       { return metadata_; }

  const Field& global_index() const { return *global_index_; }
        Field& global_index()       { return *global_index_; }

  const Field& remote_index() const { return *remote_index_; }
        Field& remote_index()       { return *remote_index_; }

  const Field& partition() const { return *partition_; }
        Field& partition()       { return *partition_; }

  const Field& ghost() const { return *ghost_; }
        Field& ghost()       { return *ghost_; }



private:

  Connectivity& add( const std::string& name, Connectivity* );

  FieldMap fields_;
  ConnectivityMap connectivities_;
  Metadata metadata_;

  // Cached shortcuts to specific fields in fields_
  Field* global_index_;
  Field* remote_index_;
  Field* partition_;
  Field* ghost_;

};

// -----------------------------------------------------------------------------------------------------

inline size_t HybridElements::size() const
{
  return size_;
}

inline size_t HybridElements::nb_types() const
{
  return element_types_.size();
}

inline const ElementType& HybridElements::element_type( size_t type_idx ) const
{
  return *element_types_[type_idx].get();
}

inline const HybridElements::Connectivity& HybridElements::node_connectivity() const
{
  return *node_connectivity_;
}

inline HybridElements::Connectivity& HybridElements::node_connectivity()
{
  return *node_connectivity_;
}

inline const HybridElements::Connectivity& HybridElements::edge_connectivity() const
{
  return *edge_connectivity_;
}

inline HybridElements::Connectivity& HybridElements::edge_connectivity()
{
  return *edge_connectivity_;
}

inline const HybridElements::Connectivity& HybridElements::cell_connectivity() const
{
  return *cell_connectivity_;
}

inline HybridElements::Connectivity& HybridElements::cell_connectivity()
{
  return *cell_connectivity_;
}


inline const Elements& HybridElements::elements( size_t type_idx ) const
{
  return *elements_[type_idx].get();
}

inline Elements& HybridElements::elements( size_t type_idx )
{
  return *elements_[type_idx].get();
}

inline size_t HybridElements::nb_nodes( size_t elem_idx ) const
{
  return node_connectivity_->rows() ? node_connectivity_->cols(elem_idx) : element_type( type_idx(elem_idx) ).nb_nodes();
}

inline size_t HybridElements::nb_edges( size_t elem_idx ) const
{
  return edge_connectivity_->rows() ? edge_connectivity_->cols(elem_idx) : element_type( type_idx(elem_idx) ).nb_edges();
}

inline size_t HybridElements::type_idx( size_t elem_idx ) const
{
  return type_idx_[elem_idx];
}

// ------------------------------------------------------------------------------------------------------

extern "C"
{
}

} // namespace mesh
} // namespace atlas

//------------------------------------------------------------------------------------------------------

#include "atlas/Parameters.h"

namespace atlas {
namespace mesh {
namespace temporary {
  
  class Convert
  {
  public:
    static HybridElements* createCells( const Mesh& );
    static HybridElements* createFaces( const Mesh& );
    static HybridElements* createElements( const Mesh&, Entity::Type );
  };
  
}

} // namespace mesh
} // namespace atlas

#endif
