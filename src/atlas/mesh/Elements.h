/*
 * (C) Copyright 1996-2016 ECMWF.
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

#ifndef atlas_Elements_H
#define atlas_Elements_H

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/array/ArrayView.h"

namespace atlas { namespace mesh { class ElementType; } }

namespace atlas {
namespace mesh {

// ------------------------------------------------------------------------------------------------------

/// @brief Describe elements of a single type
class Elements : public eckit::Owned {
public:
  typedef atlas::mesh::BlockConnectivity Connectivity;
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

  /// @brief Element to Edge connectivity table
  const Connectivity& edge_connectivity() const;
        Connectivity& edge_connectivity();

  /// @brief Element to Cell connectivity table
  const Connectivity& cell_connectivity() const;
        Connectivity& cell_connectivity();

  /// @brief Element type of these Elements
  const ElementType& element_type() const;

  /// @brief Access hybrid_elements
  /// HybridElements can contain more Elements, and holds the data.
//  const HybridElements& hybrid_elements() const;
  
  /// @brief Index of Elements in hybrid_elements
//  size_t type_idx() const;

  /// @brief Begin of elements in hybrid_elements
  size_t begin() const;

  /// @brief End of elements in hybrid_elements
  size_t end() const;
  
  
  const field::Field& field(const std::string& name) const { return hybrid_elements_->field(name); }
        field::Field& field(const std::string& name)   { return hybrid_elements_->field(name); }
  bool has_field(const std::string& name) const { return hybrid_elements_->has_field(name); }

  const field::Field& field(size_t idx) const { return hybrid_elements_->field(idx); }
        field::Field& field(size_t idx) { return hybrid_elements_->field(idx); }
  size_t nb_fields() const { return hybrid_elements_->nb_fields(); }

  const field::Field& global_index() const { return hybrid_elements_->global_index(); }
        field::Field& global_index()       { return hybrid_elements_->global_index(); }

  const field::Field& remote_index() const { return hybrid_elements_->remote_index(); }
        field::Field& remote_index()       { return hybrid_elements_->remote_index(); }

  const field::Field& partition() const { return hybrid_elements_->partition(); }
        field::Field& partition()       { return hybrid_elements_->partition(); }

  const field::Field& halo() const { return hybrid_elements_->halo(); }
        field::Field& halo()       { return hybrid_elements_->halo(); }
        
  template <typename DATATYPE, int RANK>
  array::ArrayView<DATATYPE,RANK> view( const field::Field& ) const;

  size_t add(const size_t nb_elements);

private:

  friend class HybridElements;
  void rebuild();

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

//inline size_t Elements::type_idx() const
//{
//  return type_idx_;
//}

//inline const HybridElements& Elements::hybrid_elements() const
//{
//  return *hybrid_elements_;
//}

inline const Elements::Connectivity& Elements::node_connectivity() const
{
  if( hybrid_elements_->node_connectivity().blocks() ) {
    return hybrid_elements_->node_connectivity().block(type_idx_);
  }
  else
  {
    static Elements::Connectivity dummy;
    return dummy;
  }
}

inline Elements::Connectivity& Elements::node_connectivity()
{
  if( hybrid_elements_->node_connectivity().blocks() ) {
    return hybrid_elements_->node_connectivity().block(type_idx_);
  }
  else
  {
    static Elements::Connectivity dummy;
    return dummy;
  }
}

inline const Elements::Connectivity& Elements::edge_connectivity() const
{
  if( hybrid_elements_->edge_connectivity().blocks() ) {
    return hybrid_elements_->edge_connectivity().block(type_idx_);
  }
  else
  {
    static Elements::Connectivity dummy;
    return dummy;
  }
}

inline Elements::Connectivity& Elements::edge_connectivity()
{
  if( hybrid_elements_->edge_connectivity().blocks() ) {
    return hybrid_elements_->edge_connectivity().block(type_idx_);
  }
  else
  {
    static Elements::Connectivity dummy;
    return dummy;
  }
}


inline const Elements::Connectivity& Elements::cell_connectivity() const
{
  if( hybrid_elements_->cell_connectivity().blocks() ) {
    return hybrid_elements_->cell_connectivity().block(type_idx_);
  }
  else
  {
    static Elements::Connectivity dummy;
    return dummy;
  }
}

inline Elements::Connectivity& Elements::cell_connectivity()
{
  if( hybrid_elements_->cell_connectivity().blocks() ) {
    return hybrid_elements_->cell_connectivity().block(type_idx_);
  }
  else
  {
    static Elements::Connectivity dummy;
    return dummy;
  }
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

# define field_Field field::Field
extern "C"
{
void atlas__mesh__Elements__delete(Elements* This);
size_t atlas__mesh__Elements__size(const Elements* This);
size_t atlas__mesh__Elements__begin(const Elements* This);
size_t atlas__mesh__Elements__end(const Elements* This);
BlockConnectivity* atlas__mesh__Elements__node_connectivity(Elements* This);
BlockConnectivity* atlas__mesh__Elements__edge_connectivity(Elements* This);
BlockConnectivity* atlas__mesh__Elements__cell_connectivity(Elements* This);
int atlas__mesh__Elements__has_field(const Elements* This, char* name);
int atlas__mesh__Elements__nb_fields(const Elements* This);
field_Field* atlas__mesh__Elements__field_by_idx(Elements* This, size_t idx);
field_Field* atlas__mesh__Elements__field_by_name(Elements* This, char* name);
field_Field* atlas__mesh__Elements__global_index(Elements* This);
field_Field* atlas__mesh__Elements__remote_index(Elements* This);
field_Field* atlas__mesh__Elements__partition(Elements* This);
field_Field* atlas__mesh__Elements__halo(Elements* This);
const ElementType* atlas__mesh__Elements__element_type(const Elements* This);
void atlas__mesh__Elements__add(Elements* This, size_t nb_elements);
}

#undef field_Field
//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif
