
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

PointIndex3* create_cell_centre_index( atlas::Mesh& mesh )
{
    atlas::FunctionSpace& triags = mesh.function_space( "triags" );
    ArrayView<double,2> triags_centres ( triags.field( "centre" ) );

    const std::size_t npts = triags.shape(0);

    std::vector<PointIndex3::Value> p;
    p.reserve(npts);

    for( std::size_t ip = 0; ip < npts; ++ip )
    {
        p.push_back( PointIndex3::Value(
                         PointIndex3::Point(triags_centres(ip,atlas::XX),
                                            triags_centres(ip,atlas::YY),
                                            triags_centres(ip,atlas::ZZ) ), ip ) );
    }

    PointIndex3* tree = new PointIndex3();

    tree->build(p.begin(), p.end());

    return tree;
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

