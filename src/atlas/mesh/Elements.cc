/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/runtime/ErrorHandling.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/atlas_defines.h"
#include "atlas/Field.h"

namespace atlas {
namespace mesh {

//-----------------------------------------------------------------------------

void Elements::rebuild()
{
  size_ = hybrid_elements_->elements_size_[type_idx_];
  nb_nodes_ = hybrid_elements_->element_type(type_idx_).nb_nodes();
  nb_edges_ = hybrid_elements_->element_type(type_idx_).nb_edges();
  begin_ = hybrid_elements_->elements_begin_[type_idx_];
  end_ = hybrid_elements_->elements_begin_[type_idx_+1];
}

Elements::Elements( HybridElements &elements, size_t type_idx )
  : owns_(false), hybrid_elements_(&elements), type_idx_(type_idx)
{
  rebuild();
}

Elements::Elements( ElementType* element_type, size_t nb_elements, const std::vector<idx_t> &node_connectivity )
  : owns_(true)
{
  hybrid_elements_ = new HybridElements();
  type_idx_ = hybrid_elements_->add(element_type,nb_elements,node_connectivity.data());
  rebuild();
}

Elements::Elements( ElementType* element_type, size_t nb_elements, const idx_t node_connectivity[], bool fortran_array )
  : owns_(true)
{
  hybrid_elements_ = new HybridElements();
  type_idx_ = hybrid_elements_->add(element_type,nb_elements,node_connectivity,fortran_array);
  rebuild();
}


Elements::~Elements()
{
  if( owns_ ) delete hybrid_elements_;
}

const std::string& Elements::name() const
{
  return hybrid_elements_->element_type(type_idx_).name();
}

template<> ArrayView<double,1> Elements::view( const Field& field ) const
{ 
  return ArrayView<double,1>( field.data<double>()+begin(), make_shape(size()) );
}

template<> ArrayView<float,1> Elements::view( const Field& field ) const
{
  return ArrayView<float,1>( field.data<float>()+begin(), make_shape(size()) );
}

template<> ArrayView<int,1> Elements::view( const Field& field ) const
{ 
  return ArrayView<int,1>( field.data<int>()+begin(), make_shape(size()) );
}

template<> ArrayView<long,1> Elements::view( const Field& field ) const
{ 
  return ArrayView<long,1>( field.data<long>()+begin(), make_shape(size()) );
}



template<> ArrayView<double,2> Elements::view( const Field& field ) const
{
  return ArrayView<double,2>( field.data<double>()+begin(), make_shape(size(),field.shape(1)) );
}

template<> ArrayView<float,2> Elements::view( const Field& field ) const
{
  return ArrayView<float,2>( field.data<float>()+begin(), make_shape(size(),field.shape(1)) );
}

template<> ArrayView<int,2> Elements::view( const Field& field ) const
{
  return ArrayView<int,2>( field.data<int>()+begin(), make_shape(size(),field.shape(1)) );
}

template<> ArrayView<long,2> Elements::view( const Field& field ) const
{
  return ArrayView<long,2>( field.data<long>()+begin(), make_shape(size(),field.shape(1)) );
}




size_t Elements::add(const size_t nb_elements)
{
  size_t position = size();
  hybrid_elements_->insert(end(),nb_elements);
  return position;
}

//-----------------------------------------------------------------------------

extern "C" {

void atlas__mesh__Elements__delete(Elements* This)
{
  ATLAS_ERROR_HANDLING( delete This );
}

size_t atlas__mesh__Elements__size(const Elements* This)
{
  ATLAS_ERROR_HANDLING( ASSERT( This != 0 ) );
  return This->size();
}

size_t atlas__mesh__Elements__begin(const Elements* This)
{
  ATLAS_ERROR_HANDLING( ASSERT( This != 0 ) );
  return This->begin();
}

size_t atlas__mesh__Elements__end(const Elements* This)
{
  ATLAS_ERROR_HANDLING( ASSERT( This != 0 ) );
  return This->end();
}

BlockConnectivity* atlas__mesh__Elements__node_connectivity(Elements* This)
{
  BlockConnectivity* connectivity(0);
  ATLAS_ERROR_HANDLING( connectivity = &This->node_connectivity() );
  return connectivity;
}

BlockConnectivity* atlas__mesh__Elements__edge_connectivity(Elements* This)
{
  BlockConnectivity* connectivity(0);
  ATLAS_ERROR_HANDLING( connectivity = &This->edge_connectivity() );
  return connectivity;
}

BlockConnectivity* atlas__mesh__Elements__cell_connectivity(Elements* This)
{
  BlockConnectivity* connectivity(0);
  ATLAS_ERROR_HANDLING( connectivity = &This->cell_connectivity() );
  return connectivity;
}


int atlas__mesh__Elements__has_field(const Elements* This, char* name)
{
  ATLAS_ERROR_HANDLING( ASSERT(This!=0) );
  return This->has_field(std::string(name));
}

int atlas__mesh__Elements__nb_fields(const Elements* This)
{
  ATLAS_ERROR_HANDLING( ASSERT(This!=0) );
  return This->nb_fields();
}

Field* atlas__mesh__Elements__field_by_idx(Elements* This, size_t idx)
{
  Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->field(idx);
  );
  return field;
}

Field* atlas__mesh__Elements__field_by_name(Elements* This, char* name)
{
  Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->field(std::string(name));
  );
  return field;
}

Field* atlas__mesh__Elements__global_index(Elements* This)
{
  Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->global_index();
  );
  return field;

}

Field* atlas__mesh__Elements__remote_index(Elements* This)
{
  Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->remote_index();
  );
  return field;
}

Field* atlas__mesh__Elements__partition(Elements* This)
{
  Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->partition();
  );
  return field;
}

Field* atlas__mesh__Elements__halo(Elements* This)
{
  Field* field(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    field = &This->halo();
  );
  return field;
}

const ElementType* atlas__mesh__Elements__element_type(const Elements* This)
{
  const ElementType* element_type(0);
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    element_type = &This->element_type();
  );
  return element_type;
}

void atlas__mesh__Elements__add(Elements* This, size_t nb_elements)
{
  ATLAS_ERROR_HANDLING(
    ASSERT(This!=0);
    This->add(nb_elements);
  );
}

}

//-----------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas

#undef FORTRAN_BASE
#undef TO_FORTRAN

