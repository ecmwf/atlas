/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/config/Resource.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/interpolation/method/PointIndex2.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Topology.h"

namespace atlas {
namespace interpolation {
namespace method {

ElemIndex2* create_element2D_kdtree(const Mesh& mesh, const Field& field_centres) {
    ATLAS_TRACE();
    const array::ArrayView<const double, 2> centres = array::make_view<double, 2>(field_centres);
    const array::ArrayView<const int, 1> flags      = array::make_view<int, 1>(mesh.cells().flags());
    auto include_element                            = [&](unsigned int e) {
        using util::Topology;
        return not Topology::view(flags(e)).check(Topology::INVALID);
    };

    static bool fastBuildKDTrees = eckit::Resource<bool>("$ATLAS_FAST_BUILD_KDTREES", true);

    ElemIndex2* tree      = new ElemIndex2();
    const size_t nb_cells = centres.shape(0);

    if (fastBuildKDTrees) {
        std::vector<ElemIndex2::Value> p;
        p.reserve(nb_cells);

        for (unsigned int j = 0; j < nb_cells; ++j) {
            if (include_element(j)) {
                p.emplace_back(ElemIndex2::Point(centres(j, LON), centres(j, LAT)), ElemIndex2::Payload(j));
            }
        }

        tree->build(p.begin(), p.end());
    }
    else {
        for (unsigned int j = 0; j < nb_cells; ++j) {
            if (include_element(j)) {
                tree->insert(
                    ElemIndex2::Value(ElemIndex2::Point(centres(j, LON), centres(j, LAT)), ElemIndex2::Payload(j)));
            }
        }
    }
    return tree;
}

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
