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

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildDualMesh.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"

// =======================================================

using namespace atlas::mesh::actions;

namespace atlas {
namespace numerics {
namespace fvm {

namespace {

mesh::Halo get_halo(const eckit::Parametrisation& params) {
    size_t halo_size(1);
    params.get("halo", halo_size);
    return mesh::Halo(halo_size);
}

size_t get_levels(const eckit::Parametrisation& params) {
    size_t levels(0);
    params.get("levels", levels);
    return levels;
}

double get_radius(const eckit::Parametrisation& params) {
    double radius = util::Earth::radius();
    params.get("radius", radius);
    return radius;
}

}  // namespace

Method::Method(Mesh& mesh): Method::Method(mesh, util::NoConfig()) {}

Method::Method(Mesh& mesh, const mesh::Halo& halo): Method::Method(mesh, util::Config("halo", halo.size())) {}

Method::Method(Mesh& mesh, const eckit::Configuration& params):
    mesh_(mesh),
    levels_(get_levels(params)),
    halo_(get_halo(params)),
    nodes_(mesh.nodes()),
    edges_(mesh.edges()),
    radius_(get_radius(params)) {
    setup();
}

void Method::setup() {
    ATLAS_TRACE("fvm::Method::setup ");
    util::Config config;
    config.set("halo", halo_.size());
    if (levels_) {
        config.set("levels", levels_);
    }
    node_columns_ = functionspace::NodeColumns(mesh(), config);
    edge_columns_ = functionspace::EdgeColumns(mesh(), config);

    {
        ATLAS_TRACE_SCOPE("build_median_dual_mesh") build_median_dual_mesh(mesh());
        ATLAS_TRACE_SCOPE("build_node_to_edge_connectivity") build_node_to_edge_connectivity(mesh());

        const idx_t nnodes = nodes_.size();

        auto edge_flags   = array::make_view<int, 1>(edges_.flags());
        using Topology    = mesh::Nodes::Topology;
        auto is_pole_edge = [&](size_t e) { return Topology::check(edge_flags(e), Topology::POLE); };

        // Compute sign
        {
            const mesh::Connectivity& node_edge_connectivity           = nodes_.edge_connectivity();
            const mesh::MultiBlockConnectivity& edge_node_connectivity = edges_.node_connectivity();
            if (!nodes_.has_field("node2edge_sign")) {
                nodes_.add(Field("node2edge_sign", array::make_datatype<double>(),
                                 array::make_shape(nnodes, node_edge_connectivity.maxcols())));
            }
            array::ArrayView<double, 2> node2edge_sign = array::make_view<double, 2>(nodes_.field("node2edge_sign"));

            atlas_omp_parallel_for(idx_t jnode = 0; jnode < nnodes; ++jnode) {
                for (idx_t jedge = 0; jedge < node_edge_connectivity.cols(jnode); ++jedge) {
                    idx_t iedge = node_edge_connectivity(jnode, jedge);
                    idx_t ip1   = edge_node_connectivity(iedge, 0);
                    if (jnode == ip1) {
                        node2edge_sign(jnode, jedge) = 1.;
                    }
                    else {
                        node2edge_sign(jnode, jedge) = -1.;
                        if (is_pole_edge(iedge)) {
                            node2edge_sign(jnode, jedge) = 1.;
                        }
                    }
                }
            }
        }
    }
}

// ------------------------------------------------------------------------------------------
extern "C" {
Method* atlas__numerics__fvm__Method__new(Mesh::Implementation* mesh, const eckit::Configuration* config) {
    ATLAS_ASSERT(mesh != nullptr, "Cannot access uninitialised atlas_Mesh");
    ATLAS_ASSERT(config != nullptr, "Cannot access uninitialised atlas_Config");
    Mesh m(mesh);
    return new Method(m, *config);
}

const functionspace::detail::NodeColumns* atlas__numerics__fvm__Method__functionspace_nodes(Method* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_Method");
    return dynamic_cast<const functionspace::detail::NodeColumns*>(This->node_columns().get());
}

const functionspace::detail::EdgeColumns* atlas__numerics__fvm__Method__functionspace_edges(Method* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_Method");
    return dynamic_cast<const functionspace::detail::EdgeColumns*>(This->edge_columns().get());
}
}
// ------------------------------------------------------------------------------------------

}  // namespace fvm
}  // namespace numerics
}  // namespace atlas
