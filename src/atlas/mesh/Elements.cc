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
#include "atlas/Parameters.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/Array.h"

namespace atlas {
namespace mesh {

//------------------------------------------------------------------------------------------------------

Elements::Elements(const Nodes& nodes) : nodes_(nodes) {}

Elements::~Elements() {}

void Elements::add( ElementType* element_type, size_t nb_elements )
{
  eckit::SharedPtr<ElementType> etype ( element_type );
  eckit::SharedPtr<Array> connectivity ( Array::create<int>(nb_elements,etype->nb_nodes()) );

  element_types_.push_back( etype );
  connectivity_. push_back( connectivity );
}

void Elements::add( ElementType* element_type, size_t nb_elements, const size_t connectivity[] )
{
  eckit::SharedPtr<ElementType> etype ( element_type );
  eckit::SharedPtr<Array> connectivity_array ( Array::create<int>(nb_elements,etype->nb_nodes()) );

  ArrayView<int> conn(*connectivity_array);
  size_t nb_nodes = etype->nb_nodes();
  size_t c(0);
  for( size_t e=0; e<nb_elements; ++e )
  {
    for( size_t n=0; n<nb_nodes; ++n )
    {
      conn(e,n) = connectivity[c++];
    }
  }
  element_types_.push_back( etype );
  connectivity_. push_back( connectivity_array );
}

//-----------------------------------------------------------------------------

extern "C" {

}

}  // namespace mesh
}  // namespace atlas

