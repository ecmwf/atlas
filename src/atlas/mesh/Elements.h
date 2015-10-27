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

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

namespace atlas { template<typename T> class ArrayT; }
namespace atlas { class Array; }
namespace atlas { template<typename T, int> class IndexView; }
namespace atlas { namespace mesh { class Nodes; } }
namespace atlas { namespace mesh { class ElementType; } }

namespace atlas {
namespace mesh {

class Elements;
class ElementTypeElements;

// --------------------------------------------------------------------------

class Connectivity
{
public:
  Connectivity(const Elements &elements, const int *values);
  const int* operator[](size_t elem_idx) const;
  const int* data() const { return values_; }
private:
  const int *values_;
  const Elements *elements_;
};

class ElementTypeConnectivity
{
public:
  ElementTypeConnectivity() {}
  ElementTypeConnectivity( const Elements &elements, size_t type_idx, const int *values );
  const int* operator[](size_t elem_idx) const;
  const int* data() const { return values_; }
private:
  const int *values_;
  size_t nb_nodes_;
  const Elements *elements_;
};

// -------------------------------------------------------------------------------

/**
 * \brief Elements class that describes elements in the mesh, which can are grouped
 *        per element type
 */
class Elements : public eckit::Owned {
friend class ElementTypeElements;
public: // methods

//-- Constructors

  Elements( const Nodes& nodes );
  virtual ~Elements();

//-- Accessors
  size_t size() const { return size_; }
  size_t nb_nodes(size_t elem_idx) const;
  size_t nb_edges(size_t elem_idx) const;
  const std::string& name(size_t elem_idx) const;
  size_t nb_types() const { return element_types_.size(); }
  size_t nb_elements(size_t type_idx) const { return nb_elements_[type_idx]; }
  size_t type_begin(size_t type_idx) const { return type_begin_[type_idx]; }
  size_t type_end(size_t type_idx) const { return type_end_[type_idx]; }
  size_t element_begin(size_t type_idx) const { return element_begin_[type_idx]; }
  size_t element_end(size_t type_idx) const { return element_end_[type_idx]; }
  const ElementType& element_type(size_t type_idx) const { return *element_types_[type_idx].get(); }
  const Connectivity& node_connectivity() const { return node_connectivity_access_; }
  const ElementTypeConnectivity& node_connectivity(size_t type_idx) const { return element_type_connectivity_[type_idx]; }

  // -- Modifiers

  void add( ElementType*, size_t nb_elements, const size_t node_connectivity[] );

private:
// -- Data
  const Nodes& nodes_;
  size_t size_;
  std::vector<size_t> nb_elements_;
  std::vector<size_t> type_begin_;
  std::vector<size_t> type_end_;
  std::vector< eckit::SharedPtr<ElementType> > element_types_;
  std::vector<size_t> element_begin_;
  std::vector<size_t> element_end_;
  eckit::SharedPtr< ArrayT<int> > node_connectivity_;
  eckit::SharedPtr< ArrayT<int> > nb_nodes_;
  eckit::SharedPtr< ArrayT<int> > nb_edges_;
  eckit::SharedPtr< ArrayT<int> > type_idx_;
  Connectivity node_connectivity_access_;
  std::vector<ElementTypeConnectivity> element_type_connectivity_;
};

inline const int* Connectivity::operator[](size_t elem_idx) const { return values_+elements_->element_begin(elem_idx); }

ElementTypeConnectivity::ElementTypeConnectivity(const Elements &elements, size_t type_idx, const int *values)
  : elements_(&elements),
    values_(values+elements.element_begin(elements.type_begin(type_idx))),
    nb_nodes_(elements.nb_nodes(elements.type_begin(type_idx)))
{
}

inline const int* ElementTypeConnectivity::operator[](size_t elem_idx) const {
  return values_+elem_idx*nb_nodes_;
}



class ElementTypeElements
{
public:
  ElementTypeElements(const Elements& elements, size_t type_idx);

  size_t size() const { return elements_->nb_elements(type_idx_); }
  const std::string& name() const;
  size_t nb_nodes(size_t elem_idx) const { return nb_nodes_[elem_idx]; }
  size_t nb_edges(size_t elem_idx) const { return nb_edges_[elem_idx]; }
  const ElementTypeConnectivity& node_connectivity() const { return elements_->node_connectivity(type_idx_); }

private:
  const Elements* elements_;
  size_t type_idx_;
  int* nb_nodes_;
  int* nb_edges_;
};


extern "C"
{
}

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif
