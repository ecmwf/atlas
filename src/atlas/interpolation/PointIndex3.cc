
/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/interpolation/PointIndex3.h"
#include "atlas/array/ArrayView.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace interpolation {

//------------------------------------------------------------------------------------------------------

ElemPayload make_elem_payload(size_t id, ElemPayload::ElementTypeEnum type) { return ElemPayload(id,type); }

ElemIndex3* create_element_centre_index( const atlas::mesh::Mesh& mesh )
{

    const array::ArrayView<double,2> centres ( mesh.cells().field( "centre" ) );
    const size_t ncells = mesh.cells().size();

    std::vector<ElemIndex3::Value> p;
    p.reserve(ncells);

    std::vector<ElemPayload::ElementTypeEnum> types;
    for( size_t jtype=0; jtype< mesh.cells().nb_types(); ++jtype )
      types.push_back(
            mesh.cells().element_type(jtype).name() == "Quadrilateral" ? ElemPayload::QUAD :
            mesh.cells().element_type(jtype).name() == "Triangle"      ? ElemPayload::TRIAG :
            ElemPayload::UNDEFINED
            );

    for( size_t jcell = 0; jcell < ncells; ++jcell )
    {
        p.push_back( ElemIndex3::Value(
                         ElemIndex3::Point(centres(jcell,atlas::internals::XX),
                                           centres(jcell,atlas::internals::YY),
                                           centres(jcell,atlas::internals::ZZ) ),
                         make_elem_payload(jcell,types[mesh.cells().type_idx(jcell)] ) ) );
    }


    ElemIndex3* tree = new ElemIndex3();

    tree->build(p.begin(), p.end());

    return tree;
}

//------------------------------------------------------------------------------------------------------

} // namespace interpolation
} // namespace atlas

