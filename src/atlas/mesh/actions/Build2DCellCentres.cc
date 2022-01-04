/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <limits>

#include "eckit/types/FloatCompare.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/Build2DCellCentres.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

Build2DCellCentres::Build2DCellCentres(const std::string& field_name, bool force_recompute):
    field_name_(field_name), force_recompute_(force_recompute), flatten_virtual_elements_(true) {}

Build2DCellCentres::Build2DCellCentres(eckit::Configuration& config):
    field_name_(config.getString("name", "centre")),
    force_recompute_(config.getBool("force_recompute", false)),
    flatten_virtual_elements_(config.getBool("flatten_virtual_elements", true)) {}

Field& Build2DCellCentres::operator()(Mesh& mesh) const {
    bool recompute = force_recompute_;
    if (!mesh.cells().has_field(field_name_)) {
        mesh.cells().add(Field(field_name_, array::make_datatype<double>(), array::make_shape(mesh.cells().size(), 2)));
        recompute = true;
    }
    if (recompute) {
        ATLAS_TRACE("Build2DCellCentres");
        mesh::Nodes& nodes                 = mesh.nodes();
        array::ArrayView<double, 2> coords = array::make_view<double, 2>(nodes.field("lonlat"));

        idx_t firstVirtualPoint = std::numeric_limits<idx_t>::max();
        if (nodes.metadata().has("NbRealPts")) {
            firstVirtualPoint = nodes.metadata().get<idx_t>("NbRealPts");
        }

        idx_t nb_cells = mesh.cells().size();
        std::vector<double> lons(nodes.size());
        for (idx_t e = 0; e < nodes.size(); ++e) {
            lons[e] = coords(e, LON);
            //while (lons[e] > 360.0) { lons[e] -= 360.0; }
            //while (lons[e] < 0.0) { lons[e] += 360.0; }
        }
        auto centroids = array::make_view<double, 2>(mesh.cells().field(field_name_));
        const mesh::HybridElements::Connectivity& cell_node_connectivity = mesh.cells().node_connectivity();

        for (idx_t e = 0; e < nb_cells; ++e) {
            centroids(e, LON) = 0.;
            centroids(e, LAT) = 0.;

            const idx_t nb_cell_nodes = cell_node_connectivity.cols(e);

            // check for degenerate elements (less than three unique nodes)
            // NOTE: this is not a proper check but it is very robust
            eckit::types::CompareApproximatelyEqual<double> approx(1.e-9);

            idx_t nb_equal_nodes = 0;
            for (idx_t ni = 0; ni < nb_cell_nodes - 1; ++ni) {
                idx_t i = cell_node_connectivity(e, ni);
                Point2 Pi(lons[i], coords(i, LAT));
                for (idx_t nj = ni + 1; nj < nb_cell_nodes; ++nj) {
                    idx_t j = cell_node_connectivity(e, nj);
                    Point2 Pj(lons[j], coords(j, LAT));
                    if (approx(Pi[LON], Pj[LON]) && approx(Pi[LAT], Pj[LAT])) {
                        ++nb_equal_nodes;
                    }
                }
            }

            idx_t nb_unique_nodes = idx_t(nb_cell_nodes) - nb_equal_nodes;
            if (nb_unique_nodes < 3) {
                continue;
            }

            if (flatten_virtual_elements_) {
                // calculate centroid by averaging coordinates (uses only "real" nodes)
                idx_t nb_real_nodes = 0;
                for (idx_t n = 0; n < nb_cell_nodes; ++n) {
                    const idx_t i = cell_node_connectivity(e, n);
                    if (i < firstVirtualPoint) {
                        ++nb_real_nodes;
                        centroids(e, LON) += lons[i];
                        centroids(e, LAT) += coords(i, LAT);
                    }
                }

                if (nb_real_nodes > 1) {
                    const double average_coefficient = 1. / static_cast<double>(nb_real_nodes);
                    centroids(e, LON) *= average_coefficient;
                    centroids(e, LAT) *= average_coefficient;
                }
            }
            else {
                const double average_coefficient = 1. / static_cast<double>(nb_cell_nodes);
                for (idx_t n = 0; n < nb_cell_nodes; ++n) {
                    const idx_t i = cell_node_connectivity(e, n);
                    centroids(e, LON) += lons[i] * average_coefficient;
                    centroids(e, LAT) += coords(i, LAT) * average_coefficient;
                }
            }
        }
    }
    return mesh.cells().field(field_name_);
}

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
