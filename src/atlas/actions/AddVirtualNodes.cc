/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/actions/AddVirtualNodes.h"

#include "eckit/exception/Exceptions.h"

#include "eckit/geometry/Point2.h"

#include "atlas/FunctionSpace.h"
#include "atlas/Domain.h"
#include "atlas/Mesh.h"
#include "atlas/Grid.h"

#include "atlas/grids/rgg/OctahedralRGG.h"

using eckit::geometry::LLPoint2;
using atlas::grids::rgg::OctahedralRGG;

namespace atlas {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

void AddVirtualNodes::operator()( Mesh& mesh ) const
{
    ASSERT( mesh.has_function_space("nodes") );
    ASSERT( mesh.has_grid() );

    const Grid& grid = mesh.grid();
    const Domain& domain = grid.domain();

    std::vector<LLPoint2> vPts; // virtual points

    std::vector<long> pl = OctahedralRGG::computePL(16,4); // 16 lines of latitude per hemisphere, 4 points near pole

    FunctionSpace& nodes = mesh.function_space( "nodes" );
    ArrayView<double,2> coords  ( nodes.field("xyz") );

    ArrayView<double,2> lonlat    ( nodes.field( "lonlat" ) );

}

//----------------------------------------------------------------------------------------------------------------------

} // actions
} // atlas
