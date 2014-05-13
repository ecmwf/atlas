/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"
#include "eckit/geometry/Point3.h"

#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Tesselation.h"

using namespace eckit;
using namespace eckit::geometry;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

Grid::Grid()
{
}

Grid::~Grid()
{
}

//const Mesh& Grid::mesh() const
//{
//    if( !mesh_ ) make_mesh();
//    return *mesh_;
//}

Mesh& Grid::mesh()
{
    if( !mesh_ ) make_mesh();
    return *mesh_;
}

void Grid::make_mesh()
{
    if( mesh_ ) return;

    mesh_.reset( new Mesh() );

    Mesh& mesh = *mesh_;

    const size_t npts = nPoints();

    Tesselation::create_mesh_structure( mesh, npts );

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    ASSERT(  nodes.extents()[0] == npts );

    FieldT<double>& coords  = nodes.field<double>("coordinates");
    FieldT<double>& latlon  = nodes.field<double>("latlon");
    FieldT<int>&    glb_idx = nodes.field<int>("glb_idx");

    ASSERT( npts == nodes.extents()[0] );

    const std::vector<Point>& ll = coordinates();

    for( size_t i = 0; i < npts; ++i )
    {
        glb_idx(i) = i;

        double lat = ll[i].lat();
        double lon = ll[i].lon();

        latlon(LAT,i) = lat;
        latlon(LON,i) = lon;

        eckit::geometry::latlon_to_3d( lat, lon, coords.slice(i) );
    }

    ASSERT( mesh_ );
}


//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
