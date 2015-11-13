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
#include "atlas/mesh/HybridElements.h"

namespace atlas { namespace mesh { class ElementType; } }

namespace atlas {
namespace mesh {

// ------------------------------------------------------------------------------------------------------

/// @brief Describe elements of a single type
class Elements : public eckit::Owned {
public:
  typedef atlas::BlockConnectivity Connectivity;
public:

//-- Constructors

  /// @brief Constructor that treats elements as sub-elements in HybridElements
  Elements( HybridElements &elements, size_t type_idx );

  /// @brief Constructor that internally creates a HybridElements that owns the data
  Elements( ElementType*, size_t nb_elements, const std::vector<idx_t> &node_connectivity );

  /// @brief Constructor that internally creates a HybridElements that owns the data
  Elements( ElementType*, size_t nb_elements, const idx_t node_connectivity[], bool fortran_array=false );

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

  /// @brief Begin of elements in hybrid_elements
  size_t begin() const;

  /// @brief End of elements in hybrid_elements
  size_t end() const;

private:
  bool owns_;
  HybridElements* hybrid_elements_;
  size_t size_;
  size_t begin_;
  size_t end_;
  size_t type_idx_;
  size_t nb_nodes_;
  size_t nb_edges_;
};

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

inline size_t Elements::begin() const
{
  return begin_;
}

inline size_t Elements::end() const
{
  return end_;
}

// ------------------------------------------------------------------------------------------------------

extern "C"
{
}

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif
