/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/mesh/actions/AddVirtualNodes.h"

#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "atlas/field/Field.h"
#include "atlas/grid/domains/RectangularDomain.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/gaussian/OctahedralGaussian.h"
#include "atlas/internals/Parameters.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"


namespace atlas {
namespace mesh {
namespace actions {


void AddVirtualNodes::operator()(const atlas::grid::Grid& grid, atlas::mesh::Mesh& mesh) const {
    using eckit::geometry::LLPoint2;
    const grid::Domain& dom = grid.domain();

    if (dom.isGlobal()) return; // don't add virtual points to global domains

    const grid::Grid& octa = atlas::grid::gaussian::OctahedralGaussian(16);

    std::vector<LLPoint2> allPts;
    octa.lonlat(allPts);

    std::vector<LLPoint2> vPts; // virtual points

    // loop over the point and keep the ones that *don't* fall in the domain
    for(size_t i = 0; i < allPts.size(); ++i) {
        const LLPoint2& p = allPts[i];
        if( !dom.contains(p.lon(),p.lat()) )
            vPts.push_back(p);
    }

    mesh::Nodes& nodes = mesh.nodes();

    const size_t nb_real_pts = nodes.size();
    const size_t nb_virtual_pts = vPts.size();

    size_t new_size = nodes.size() + vPts.size();

    nodes.resize(new_size); // resizes the fields

    const size_t nb_total_pts = nodes.size();

    ASSERT( nb_total_pts == nb_real_pts + nb_virtual_pts );

    nodes.metadata().set<size_t>("NbRealPts",nb_real_pts);
    nodes.metadata().set<size_t>("NbVirtualPts",nb_virtual_pts);

    array::ArrayView<double,2> coords ( nodes.field("xyz") );
    array::ArrayView<double,2> lonlat ( nodes.lonlat() );
    array::ArrayView<gidx_t,1> gidx   ( nodes.global_index() );

    for(size_t i = 0; i < nb_virtual_pts; ++i) {
        const size_t n = nb_real_pts + i;
        lonlat(n,internals::LON) = vPts[i].lon();
        lonlat(n,internals::LAT) = vPts[i].lat();
        eckit::geometry::lonlat_to_3d(lonlat[n].data(),coords[n].data());
        gidx(n) = n+1;
    }
}


} // namespace actions
} // namespace mesh
} // namespace atlas
