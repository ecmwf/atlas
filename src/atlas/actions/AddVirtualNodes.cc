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
#include "eckit/types/Types.h"
#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"

#include "atlas/FunctionSpace.h"
#include "atlas/Domain.h"
#include "atlas/Mesh.h"
#include "atlas/Nodes.h"
#include "atlas/Grid.h"

#include "atlas/grids/rgg/OctahedralRGG.h"

using eckit::Log;
using eckit::geometry::LLPoint2;
using atlas::grids::rgg::OctahedralRGG;
using eckit::operator<<;

namespace atlas {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

void AddVirtualNodes::operator()( Mesh& mesh ) const
{
    ASSERT( mesh.has_function_space("nodes") );
    ASSERT( mesh.has_grid() );

    const Grid& grid = mesh.grid();
    const Domain& domain = grid.domain();

    if( domain.global() ) return; // don't add virtual points to global domains

    const Grid& octa = OctahedralRGG(16,4);

    std::vector<LLPoint2> allPts;
    octa.lonlat(allPts);

    std::vector<LLPoint2> vPts; // virtual points

    // loop over the point and keep the ones that *don't* fall in the domain
    for(size_t i = 0; i < allPts.size(); ++i)
    {
        const LLPoint2& p = allPts[i];
        if( !domain.contains(p.lon(),p.lat()) )
            vPts.push_back(p);
    }

    Nodes& nodes = mesh.nodes();

    const size_t nb_real_pts = nodes.size();
    const size_t nb_virtual_pts = vPts.size();

    size_t new_size = nodes.size() + vPts.size();

    nodes.resize(new_size); // resizes the fields

    const size_t nb_total_pts = nodes.size();

    ASSERT( nb_total_pts == nb_real_pts + nb_virtual_pts );

    nodes.metadata().set<size_t>("NbRealPts",nb_real_pts);
    nodes.metadata().set<size_t>("NbVirtualPts",nb_virtual_pts);

    ArrayView<double,2> coords ( nodes.field("xyz") );
    ArrayView<double,2> lonlat ( nodes.lonlat() );
    ArrayView<gidx_t,1> gidx   ( nodes.global_index() );

    for(size_t i = 0; i < nb_virtual_pts; ++i)
    {
        const size_t n = nb_real_pts + i;
        lonlat(n,LON) = vPts[i].lon();
        lonlat(n,LAT) = vPts[i].lat();
        eckit::geometry::lonlat_to_3d(lonlat[n].data(),coords[n].data());
        gidx(n) = n+1;
    }
}

//----------------------------------------------------------------------------------------------------------------------

} // actions
} // atlas
