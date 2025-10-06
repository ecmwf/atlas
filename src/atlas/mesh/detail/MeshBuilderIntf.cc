/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/detail/MeshBuilderIntf.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace mesh {

//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

TriangularMeshBuilder* atlas__TriangularMeshBuilder__new() {
    return new TriangularMeshBuilder();
}

void atlas__TriangularMeshBuilder__delete(TriangularMeshBuilder* This) {
    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_TriangularMeshBuilder");
    delete This;
}

Mesh::Implementation* atlas__TriangularMeshBuilder__operator(TriangularMeshBuilder* This,
        size_t nb_nodes, const gidx_t node_global_index[],
        const double x[], const double y[], size_t xstride, size_t ystride,
        const double lon[], const double lat[], size_t lonstride, size_t latstride,
        size_t nb_triags, const gidx_t triangle_global_index[], const gidx_t triangle_nodes_global_index[],
        gidx_t global_index_base) {

    ATLAS_ASSERT(This != nullptr, "Cannot access uninitialisd atlas_TriangularMeshBuilder");

    Mesh::Implementation* m;
    {
        Mesh mesh = This->operator()(nb_nodes, node_global_index,
                                     x, y, xstride, ystride,
                                     lon, lat, lonstride, latstride,
                                     nb_triags, triangle_global_index, triangle_nodes_global_index,
                                     global_index_base);
        mesh.get()->attach();
        m = mesh.get();
    }
    m->detach();
    return m;
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
