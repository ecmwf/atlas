/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grid/PointIndex3.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

PointIndex3* create_cell_centre_index( atlas::Mesh& mesh )
{
    atlas::FunctionSpace& triags = mesh.function_space( "triags" );
    atlas::FieldT<double>& triags_centres = triags.field<double>( "centre" );

    const size_t npts = triags.extents()[0];

    std::vector<typename PointIndex3::Value> p;
    p.reserve(npts);

    for( size_t ip = 0; ip < npts; ++ip )
    {
        p.push_back( typename PointIndex3::Value(
                         typename PointIndex3::Point( triags_centres(atlas::XX,ip),
                                                      triags_centres(atlas::YY,ip),
                                                      triags_centres(atlas::ZZ,ip) ), ip ) );
    }

    PointIndex3* tree = new PointIndex3();

    tree->build(p.begin(), p.end());

    return tree;
}

//------------------------------------------------------------------------------------------------------

} // namespace atlas

