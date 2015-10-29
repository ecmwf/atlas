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
/// @date October 2015

#ifndef atlas_mesh_Elements_H
#define atlas_mesh_Elements_H

#include "atlas/atlas_config.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/util/Debug.h"

namespace atlas { template<typename T> class ArrayT; }
namespace atlas { namespace mesh { class ElementType; } }

namespace atlas {
namespace mesh {

// Classes defined in this file:
class HybridElements;
class Elements;
class HybridConnectivity;
class BlockConnectivity;

// --------------------------------------------------------------------------

#ifdef ATLAS_HAVE_FORTRAN
#define FROM_FORTRAN -1
#define TO_FORTRAN +1
#else
#define FROM_FORTRAN
#define TO_FORTRAN
#endif

class HybridConnectivity : public eckit::Owned
{
public:
  HybridConnectivity( idx_t *values, size_t rows, size_t *row_offset, size_t blocks, size_t *block_offset );
  idx_t operator()(size_t row, size_t col) const;
  void set(size_t row, const idx_t column_values[]);
  size_t rows() const { return rows_; }
  size_t cols(size_t row) const { return row_offset_[row+1]-row_offset_[row]; }
  size_t blocks() const { return blocks_; }

  HybridConnectivity();
  ~HybridConnectivity();

  void add( size_t rows, size_t cols, const idx_t values[], bool fortran_array=false );
  void add( const BlockConnectivity& );

  const BlockConnectivity& block_connectivity(size_t block_idx) const { return *block_connectivity_[block_idx].get(); }
        BlockConnectivity& block_connectivity(size_t block_idx)       { return *block_connectivity_[block_idx].get(); }

private:
  void regenerate_block_connectivity();

private:
  bool owns_;
  std::vector<idx_t>  owned_values_;
  std::vector<size_t> owned_row_offset_;
  std::vector<size_t> owned_block_offset_;

  idx_t  *values_;
  size_t rows_;
  size_t *row_offset_;
  size_t blocks_;
  size_t *block_offset_;

  std::vector< eckit::SharedPtr<BlockConnectivity> > block_connectivity_;
};

class BlockConnectivity : public eckit::Owned
{
  friend class HybridConnectivity;
public:
  BlockConnectivity();
  BlockConnectivity( size_t rows, size_t cols, idx_t *values );
  idx_t operator()( size_t row, size_t col ) const;
  void set( size_t row, const idx_t column_values[] );
  void add( size_t rows, size_t cols, const idx_t *values, bool fortran_array=false );
  size_t rows() const { return rows_; }
  size_t cols() const { return cols_; }
private:
  bool owns_;
  std::vector<idx_t> owned_values_;

  size_t rows_;
  size_t cols_;
  idx_t *values_;

};

// -------------------------------------------------------------------------------

class Elements
{
public:
  typedef atlas::mesh::BlockConnectivity Connectivity;
public:
  Elements();

  // Constructor that treats elements as sub-elements in HybridElements
  Elements(HybridElements &elements, size_t type_idx);

  // Constructor that internally creates a HybridElements
  Elements(ElementType*, size_t nb_elements, const std::vector<idx_t> &node_connectivity );
  Elements(ElementType*, size_t nb_elements, const idx_t node_connectivity[], bool fortran_array=false );

  virtual ~Elements();

  size_t size() const;
  const std::string& name() const;
  size_t nb_nodes() const;
  size_t nb_edges() const;
  const Connectivity& node_connectivity() const;
  void set_node_connectivity( size_t elem_idx, const idx_t node_connectivity[] );
  const ElementType& element_type() const;

  const HybridElements& hybrid_elements() const { return *hybrid_elements_; }

private:
  HybridElements* hybrid_elements_;
  size_t type_idx_;
  size_t nb_nodes_;
  size_t nb_edges_;
  bool owns_elements_;
};

// -------------------------------------------------------------------------------

/**
 * \brief HybridElements class that describes elements in the mesh, which can are grouped
 *        per element type
 */
class HybridElements : public eckit::Owned {
friend class Elements;
public:
  typedef HybridConnectivity Connectivity;

public: // methods

//-- Constructors

  HybridElements();
  virtual ~HybridElements();

//-- Accessors

  size_t size() const { return size_; }
  size_t nb_nodes(size_t hybrid_elem_idx) const;
  size_t nb_edges(size_t hybrid_elem_idx) const;
  const std::string& name(size_t hybrid_elem_idx) const;
  size_t nb_types() const;
  const ElementType& element_type(size_t type_idx) const;
  const HybridElements::Connectivity& node_connectivity() const;
  const Elements::Connectivity& node_connectivity(size_t type_idx) const;
  const Elements& elements(size_t type_idx) const;
        Elements& elements(size_t type_idx);

  // Advanced api. to be seen if needed
  //  size_t nb_elements(size_t type_idx) const { return elements_size_[type_idx]; }
  //  size_t type_begin(size_t type_idx) const { return elements_begin_[type_idx]; }
  //  size_t type_end(size_t type_idx) const { return elements_end_[type_idx]; }
  //  size_t element_nodes_begin(size_t type_idx) const { return nodes_begin_[type_idx]; }
  //  size_t element_end(size_t type_idx) const { return nodes_end_[type_idx]; }

// -- Modifiers

  size_t add( const ElementType*, size_t nb_elements, const std::vector<idx_t> &node_connectivity );
  size_t add( const ElementType*, size_t nb_elements, const idx_t node_connectivity[] );
  size_t add( const ElementType*, size_t nb_elements, const idx_t node_connectivity[], bool fortran_array );
  size_t add( const Elements& );

  void set_node_connectivity( size_t elem_idx, const idx_t node_connectivity[] );
  void set_node_connectivity( size_t type_idx, size_t elem_idx, const idx_t node_connectivity[] );

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
  eckit::SharedPtr< ArrayT<idx_t> > node_connectivity_array_;

// -- Accessor helpers
  eckit::SharedPtr< HybridElements::Connectivity > node_connectivity_;
  std::vector< eckit::SharedPtr<Elements::Connectivity> > node_block_connectivity_;
  std::vector<Elements> elements_;
};

// -----------------------------------------------------------------------------------------------------

inline idx_t HybridConnectivity::operator()(size_t row, size_t col) const
{
  return (values_+row_offset_[row])[col] FROM_FORTRAN;
}

inline void HybridConnectivity::set(size_t row, const idx_t column_values[]) {
  idx_t *col = values_+row_offset_[row];
  size_t N = row_offset_[row+1]-row_offset_[row];
  for( size_t n=0; n<N; ++n ) {
    col[n] = column_values[n] TO_FORTRAN;
  }
}

// -----------------------------------------------------------------------------------------------------

BlockConnectivity::BlockConnectivity( size_t rows, size_t cols, idx_t *values )
  : rows_(rows),
    cols_(cols),
    values_(values)
{
}

inline idx_t BlockConnectivity::operator()(size_t row, size_t col) const {
  return (values_+row*cols_)[col] FROM_FORTRAN;
}

inline void BlockConnectivity::set(size_t row, const idx_t column_values[]) {
  DEBUG("set");
  idx_t *col = values_+row*cols_;
  DEBUG_VAR(row);
  DEBUG_VAR(row*cols_);
  DEBUG_VAR(rows_*cols_);
  for( size_t n=0; n<cols_; ++n ) {
    DEBUG_VAR(n);
    col[n] = column_values[n] TO_FORTRAN;
  }
}

// ------------------------------------------------------------------------------------------------------

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

inline const Elements::Connectivity& HybridElements::node_connectivity(size_t type_idx) const
{
  return node_connectivity_.get()->block_connectivity(type_idx);
}

inline const Elements& HybridElements::elements(size_t type_idx) const
{
  return elements_[type_idx];
}

inline Elements& HybridElements::elements(size_t type_idx)
{
  return elements_[type_idx];
}

// ------------------------------------------------------------------------------------------------------

inline size_t Elements::size() const
{ return hybrid_elements_->elements_size_[type_idx_]; }

inline size_t Elements::nb_nodes() const
{ return nb_nodes_; }

inline size_t Elements::nb_edges() const
{ return nb_edges_; }

inline const Elements::Connectivity& Elements::node_connectivity() const
{
  return hybrid_elements_->node_connectivity(type_idx_);
}
inline const ElementType& Elements::element_type() const
{
  return hybrid_elements_->element_type(type_idx_);
}

// ------------------------------------------------------------------------------------------------------

extern "C"
{
}

#undef FROM_FORTRAN
#undef TO_FORTRAN

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif
