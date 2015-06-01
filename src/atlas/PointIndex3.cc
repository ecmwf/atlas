
/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/util/ArrayView.h"
#include "atlas/PointIndex3.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

ElemPayload make_elem_payload(size_t id, char t) { return ElemPayload(id,t); }

ElemIndex3* create_element_centre_index( const atlas::Mesh& mesh )
{
    atlas::FunctionSpace& triags = mesh.function_space( "triags" );
    atlas::FunctionSpace& quads  = mesh.function_space( "quads" );

    ArrayView<double,2> triags_centres ( triags.field( "centre" ) );
    ArrayView<double,2> quads_centres  ( quads.field( "centre" ) );

    const size_t ntriags = triags.shape(0);
    const size_t nquads  = quads.shape(0);

    std::vector<ElemIndex3::Value> p;
    p.reserve(ntriags+nquads);

    for( size_t ip = 0; ip < ntriags; ++ip )
    {
        p.push_back( ElemIndex3::Value(
                         ElemIndex3::Point(triags_centres(ip,atlas::XX),
                                           triags_centres(ip,atlas::YY),
                                           triags_centres(ip,atlas::ZZ) ), make_elem_payload(ip,'t') ) );
    }

    for( size_t ip = 0; ip < nquads; ++ip )
    {
        p.push_back( ElemIndex3::Value(
                         ElemIndex3::Point(quads_centres(ip,atlas::XX),
                                           quads_centres(ip,atlas::YY),
                                           quads_centres(ip,atlas::ZZ) ), make_elem_payload(ip,'q') ) );
    }

    ElemIndex3* tree = new ElemIndex3();

    tree->build(p.begin(), p.end());

    return tree;
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

