/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include "eckit/log/BigNum.h"
#include "atlas/util/QhullSphericalTriangulation.h"
#include "atlas/util/CGALSphericalTriangulation.h"
#include "atlas/library/config.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Grid.h"
#include "atlas/interpolation/method/PointSet.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildConvexHull3D.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

using atlas::interpolation::method::PointSet;

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

BuildConvexHull3D::BuildConvexHull3D(const eckit::Parametrisation& config) {
    config.get("remove_duplicate_points", remove_duplicate_points_ = true);
}

//----------------------------------------------------------------------------------------------------------------------

void BuildConvexHull3D::operator()(Mesh& mesh) const {
    // don't tesselate meshes already with triags or quads
    if (mesh.cells().size()) {
        return;
    }

    std::string default_backend = (ATLAS_HAVE_CGAL ? "cgal" : "qhull");
    std::string backend = eckit::Resource<std::string>("$ATLAS_DELAUNAY_BACKEND",default_backend);

    ATLAS_TRACE("BuildConvexHull3D [" + backend + "]");

    std::vector<size_t> local_index;
    if (remove_duplicate_points_) {
        PointSet points(mesh);
        points.list_unique_points(local_index);
    }
    
    auto add_triangles = [&]( const auto& triangles ) {
        const idx_t nb_triags = triangles.size();
        mesh.cells().add(mesh::ElementType::create("Triangle"), nb_triags);
        mesh::HybridElements::Connectivity& triag_nodes = mesh.cells().node_connectivity();
        auto triag_gidx = array::make_view<gidx_t, 1>(mesh.cells().global_index());
        auto triag_part = array::make_view<int, 1>(mesh.cells().partition());

        Log::debug() << "Inserting triags (" << eckit::BigNum(nb_triags) << ")" << std::endl;

        for (idx_t tidx = 0; tidx<nb_triags; ++tidx) {
            auto& t = triangles[tidx];
            std::array<idx_t,3> idx{t[0],t[1],t[2]};
            if( local_index.size() ) {
                idx[0] = local_index[idx[0]];
                idx[1] = local_index[idx[1]];
                idx[2] = local_index[idx[2]];
            }
            triag_nodes.set(tidx, idx.data());
            triag_gidx(tidx) = tidx + 1;
            triag_part(tidx) = 0;
        }
    };


    if( local_index.size() == mesh.nodes().size() or local_index.empty() ) {
        local_index.clear();
        auto lonlat = array::make_view<double,2>(mesh.nodes().lonlat());
        if (backend == "qhull") {
            ATLAS_TRACE("qhull");
            auto triangles = util::QhullSphericalTriangulation{static_cast<size_t>(lonlat.shape(0)),lonlat.data()}.triangles();
            add_triangles(triangles);
        }
        else if (backend == "cgal") {
            ATLAS_TRACE("cgal");
            auto triangles = util::CGALSphericalTriangulation{static_cast<size_t>(lonlat.shape(0)),lonlat.data()}.triangles();
            add_triangles(triangles);
        }
        else {
            ATLAS_THROW_EXCEPTION("backend " << backend << " not supported");
        }
    }
    else {
        auto lonlat_view = array::make_view<double,2>(mesh.nodes().lonlat());

        std::vector<PointLonLat> lonlat(local_index.size());
        size_t jnode = 0;
        for( auto& ip: local_index ) {
            lonlat[jnode] = {lonlat_view(ip,0),lonlat_view(ip,1)};
            ++jnode;
        }
        if (backend == "qhull") {
            ATLAS_TRACE("qhull");
            auto triangles = util::QhullSphericalTriangulation{lonlat}.triangles();
            add_triangles(triangles);
        }
        else if (backend == "cgal") {
            ATLAS_TRACE("cgal");
            auto triangles = util::CGALSphericalTriangulation{lonlat}.triangles();
            add_triangles(triangles);
        }
        else {
            ATLAS_THROW_EXCEPTION("backend " << backend << " not supported");
        }
    }

}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
