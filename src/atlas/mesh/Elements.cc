/*
 * (C) Copyright 1996-2015 ECMWF.
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

}

//-----------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas

#undef FORTRAN_BASE
#undef TO_FORTRAN

