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

namespace atlas { template<typename T> class ArrayT; }
namespace atlas { class Array; }
namespace atlas { template<typename T, int> class IndexView; }
namespace atlas { namespace mesh { class Nodes; } }
namespace atlas { namespace mesh { class ElementType; } }

namespace atlas {
namespace mesh {

// Classes defined in this file:
class HybridElements;
class Elements;
class HybridConnectivity;
class Connectivity;

// --------------------------------------------------------------------------

#ifdef ATLAS_HAVE_FORTRAN
#define FROM_FORTRAN -1
#define TO_FORTRAN +1
#else
#define FROM_FORTRAN
#define TO_FORTRAN
#endif

class HybridConnectivity
{
public:
  HybridConnectivity(idx_t *values, size_t *offset);
  idx_t operator()(size_t row, size_t col) const;
private:
  idx_t  *values_;
  size_t *offset_;
};

class Connectivity
{
public:
  Connectivity() {}
  Connectivity( idx_t *values, size_t stride );
  idx_t operator()(size_t row, size_t col) const;
  void set(size_t row, const idx_t column_values[]);
private:
  idx_t *values_;
  size_t stride_;
};

// -------------------------------------------------------------------------------

class Elements
{
public:
  typedef atlas::mesh::Connectivity Connectivity;
public:
  Elements();

  // Constructor that treats elements as sub-elements in HybridElements
  Elements(HybridElements &elements, size_t type_idx);

  // Constructor that internally creates a HybridElements
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
  size_t nb_nodes(size_t elem_idx) const;
  size_t nb_edges(size_t elem_idx) const;
  const std::string& name(size_t elem_idx) const;
  size_t nb_types() const { return element_types_.size(); }
  const ElementType& element_type(size_t type_idx) const { return *element_types_[type_idx].get(); }
  const HybridElements::Connectivity& node_connectivity() const { return node_connectivity_access_; }
  const Elements::Connectivity& node_connectivity(size_t type_idx) const { return element_type_connectivity_[type_idx]; }
  const Elements& elements(size_t type_idx) const { return elements_[type_idx]; }
        Elements& elements(size_t type_idx)       { return elements_[type_idx]; }

  // Advanced api. to be seen if needed
  //  size_t nb_elements(size_t type_idx) const { return nb_elements_[type_idx]; }
  //  size_t type_begin(size_t type_idx) const { return type_begin_[type_idx]; }
  //  size_t type_end(size_t type_idx) const { return type_end_[type_idx]; }
  //  size_t element_begin(size_t type_idx) const { return element_begin_[type_idx]; }
  //  size_t element_end(size_t type_idx) const { return element_end_[type_idx]; }

// -- Modifiers

  size_t add( const ElementType*, size_t nb_elements, const idx_t node_connectivity[] );
  size_t add( const ElementType*, size_t nb_elements, const idx_t node_connectivity[], bool fortran_array );
  size_t add( const Elements& );

  void set_node_connectivity( size_t elem_idx, const idx_t node_connectivity[] );
  void set_node_connectivity( size_t type_idx, size_t elem_idx, const idx_t node_connectivity[] );

private:

// -- Data
  size_t size_;                          //!< total number of elements

// -- Data: one value per type
  std::vector<size_t> nb_elements_;
  std::vector<size_t> type_begin_;
  std::vector<size_t> type_end_;
  std::vector< eckit::SharedPtr<const ElementType> > element_types_;

// -- Data: one value per element
  std::vector<size_t> element_begin_;
  std::vector<size_t> element_end_;
  std::vector<size_t> nb_nodes_;
  std::vector<size_t> nb_edges_;
  std::vector<size_t> type_idx_;

// -- Data: one value per node per element
  eckit::SharedPtr< ArrayT<idx_t> > node_connectivity_;

// -- Accessor helpers
  HybridElements::Connectivity node_connectivity_access_;
  std::vector<Elements::Connectivity> element_type_connectivity_;
  std::vector<Elements> elements_;
};

// -----------------------------------------------------------------------------------------------------

inline idx_t HybridConnectivity::operator()(size_t row, size_t col) const
{
  return (values_+offset_[row])[col] FROM_FORTRAN;
}

// -----------------------------------------------------------------------------------------------------

Connectivity::Connectivity(idx_t *values, size_t stride)
  : values_(values),
    stride_(stride)
{
}

inline idx_t Connectivity::operator()(size_t row, size_t col) const {
  return (values_+row*stride_)[col] FROM_FORTRAN;
}

inline void Connectivity::set(size_t row, const idx_t column_values[]) {
  idx_t *col = values_+row*stride_;
  for( size_t n=0; n<stride_; ++n ) {
    col[n] = column_values[n] TO_FORTRAN;
  }
}



inline size_t Elements::size() const { return hybrid_elements_->nb_elements_[type_idx_]; }
inline size_t Elements::nb_nodes() const { return nb_nodes_; }
inline size_t Elements::nb_edges() const { return nb_edges_; }
inline const Elements::Connectivity& Elements::node_connectivity() const { return hybrid_elements_->node_connectivity(type_idx_); }
inline const ElementType& Elements::element_type() const { return hybrid_elements_->element_type(type_idx_); }



extern "C"
{
}

#undef FROM_FORTRAN
#undef TO_FORTRAN

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif
