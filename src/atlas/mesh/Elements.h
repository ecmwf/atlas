/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @file Elements.h
/// @author Willem Deconinck
/// @date October 2015
///
/// This file describes the Elements class for a Mesh.

#ifndef atlas_mesh_Elements_H
#define atlas_mesh_Elements_H

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/Connectivity.h"

namespace atlas { namespace mesh { class ElementType; } }

namespace atlas {
namespace mesh {

// Classes defined in this file:
class HybridElements;
class Elements;

// -------------------------------------------------------------------------------

/// @brief Describe elements of a single type
class Elements : public eckit::Owned {
public:
  typedef atlas::BlockConnectivity Connectivity;
public:

//-- Constructors

  /// @brief Constructor that treats elements as sub-elements in HybridElements
  Elements(HybridElements &elements, size_t type_idx);

  /// @brief Constructor that internally creates a HybridElements that owns the data
  Elements(ElementType*, size_t nb_elements, const std::vector<idx_t> &node_connectivity );

  /// @brief Constructor that internally creates a HybridElements that owns the data
  Elements(ElementType*, size_t nb_elements, const idx_t node_connectivity[], bool fortran_array=false );

  /// @brief Destructor
  virtual ~Elements();

//-- Accessors

  /// @brief Number of elements
  size_t size() const;

  /// @brief Name of this element type
  const std::string& name() const;

  /// @brief Number of nodes for each element type
  size_t nb_nodes() const;

  /// @brief Number of edges for each element type
  size_t nb_edges() const;

  /// @brief Element to Node connectivity table
  const Connectivity& node_connectivity() const;
        Connectivity& node_connectivity();

  /// @brief Element type of these Elements
  const ElementType& element_type() const;

  /// @brief Access hybrid_elements
  /// HybridElements can contain more Elements, and holds the data.
  const HybridElements& hybrid_elements() const { return *hybrid_elements_; }

private:
  bool owns_;
  HybridElements* hybrid_elements_;
  size_t size_;
  size_t type_idx_;
  size_t nb_nodes_;
  size_t nb_edges_;
};

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
  size_t nb_nodes(size_t elem_idx) const;

  /// @brief Number of edges for given element
  size_t nb_edges(size_t elem_idx) const;

  /// @brief Element type index for given element
  size_t type_idx(size_t elem_idx) const;

  /// @brief Element type name for given element
  const std::string& name(size_t elem_idx) const;

  /// @brief Element to Node connectivity table
  const HybridElements::Connectivity& node_connectivity() const;
        HybridElements::Connectivity& node_connectivity();

  /// @brief Number of types present in HybridElements
  size_t nb_types() const;

  /// @brief The element_type description for given type
  const ElementType& element_type(size_t type_idx) const;

  /// @brief Sub-elements convenience class for given type
  /// This allows optimized access to connectivities and loops.
  const Elements& elements(size_t type_idx) const;
        Elements& elements(size_t type_idx);

// -- Modifiers

  /// @brief Add a new element type with given number of elements and node-connectivity
  size_t add( const ElementType*, size_t nb_elements, const std::vector<idx_t> &node_connectivity );

  /// @brief Add a new element type with given number of elements and node-connectivity
  size_t add( const ElementType*, size_t nb_elements, const idx_t node_connectivity[] );

  /// @brief Add a new element type with given number of elements and node-connectivity
  size_t add( const ElementType*, size_t nb_elements, const idx_t node_connectivity[], bool fortran_array );

  /// @brief Add a new element type from existing Elements.
  /// Data will be copied.
  size_t add( const Elements& );

private:

// -- Data
  size_t size_;                          //!< total number of elements

// -- Data: one value per type
  std::vector<size_t> elements_size_;
  std::vector<size_t> elements_begin_;
  std::vector< eckit::SharedPtr<const ElementType> > element_types_;

// -- Data: one value per element
  std::vector<size_t> nodes_begin_;
  std::vector<size_t> nb_nodes_;
  std::vector<size_t> nb_edges_;
  std::vector<size_t> type_idx_;

// -- Data: one value per node per element
  std::vector<idx_t> node_connectivity_array_;

// -- Accessor helpers
  eckit::SharedPtr< HybridElements::Connectivity > node_connectivity_;
  std::vector< eckit::SharedPtr<Elements> > elements_;
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

inline const ElementType& HybridElements::element_type(size_t type_idx) const
{
  return *element_types_[type_idx].get();
}

inline const HybridElements::Connectivity& HybridElements::node_connectivity() const
{
  return *node_connectivity_.get();
}

inline HybridElements::Connectivity& HybridElements::node_connectivity()
{
  return *node_connectivity_.get();
}

inline const Elements& HybridElements::elements(size_t type_idx) const
{
  return *elements_[type_idx].get();
}

inline Elements& HybridElements::elements(size_t type_idx)
{
  return *elements_[type_idx].get();
}

inline size_t HybridElements::nb_nodes(size_t elem_idx) const
{
  return nb_nodes_[elem_idx];
}

inline size_t HybridElements::nb_edges(size_t elem_idx) const
{
  return nb_edges_[elem_idx];
}

inline size_t HybridElements::type_idx(size_t elem_idx) const
{
  return type_idx_[elem_idx];
}


// ------------------------------------------------------------------------------------------------------

inline size_t Elements::size() const
{
  return size_;
}

inline size_t Elements::nb_nodes() const
{
  return nb_nodes_;
}

inline size_t Elements::nb_edges() const
{
  return nb_edges_;
}

inline const Elements::Connectivity& Elements::node_connectivity() const
{
  return hybrid_elements_->node_connectivity().block(type_idx_);
}

inline Elements::Connectivity& Elements::node_connectivity()
{
  return hybrid_elements_->node_connectivity().block(type_idx_);
}

inline const ElementType& Elements::element_type() const
{
  return hybrid_elements_->element_type(type_idx_);
}

// ------------------------------------------------------------------------------------------------------

extern "C"
{
}

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif
