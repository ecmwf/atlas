/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <limits>
#include "eckit/types/FloatCompare.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/field/Field.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/array/MakeView.h"

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

BuildCellCentres::BuildCellCentres( const std::string& field_name ) :
    field_name_(field_name) {
}

Field& BuildCellCentres::operator()( Mesh& mesh ) const {
    if (!mesh.cells().has_field(field_name_)) {

        mesh::Nodes& nodes = mesh.nodes();
        array::ArrayView<double, 2> coords = array::make_view<double, 2>(nodes.field("xyz"));

        size_t firstVirtualPoint = std::numeric_limits<size_t>::max();
        if (nodes.metadata().has("NbRealPts")) {
            firstVirtualPoint = nodes.metadata().get<size_t>("NbRealPts");
        }

        size_t nb_cells = mesh.cells().size();
        array::ArrayView<double, 2> centroids = array::make_view<double, 2>(
                    mesh.cells().add(
                        Field(field_name_, array::make_datatype<double>(), array::make_shape(nb_cells, 3)) ));
        const mesh::HybridElements::Connectivity& cell_node_connectivity = mesh.cells().node_connectivity();

        for (size_t e = 0; e < nb_cells; ++e) {
            centroids(e, XX) = 0.;
            centroids(e, YY) = 0.;
            centroids(e, ZZ) = 0.;

            const size_t nb_cell_nodes = cell_node_connectivity.cols(e);

            // check for degenerate elements (less than three unique nodes)
            // NOTE: this is not a proper check but it is very robust
            eckit::types::CompareApproximatelyEqual<double> approx(1.e-9);

            int nb_equal_nodes = 0;
            for (size_t ni = 0; ni < nb_cell_nodes - 1; ++ni) {
                idx_t i = cell_node_connectivity(e, ni);
                Point3 Pi(coords(i, XX), coords(i, YY), coords(i, ZZ));
                for (size_t nj = ni + 1; nj < nb_cell_nodes; ++nj) {
                    idx_t j = cell_node_connectivity(e, nj);
                    Point3 Pj(coords(j, XX), coords(j, YY), coords(j, ZZ));
                    if (approx(Pi[XX], Pj[XX]) &&
                        approx(Pi[YY], Pj[YY]) &&
                        approx(Pi[ZZ], Pj[ZZ])) {
                        ++nb_equal_nodes;
                    }
                }
            }

            int nb_unique_nodes = int(nb_cell_nodes) - nb_equal_nodes;
            if (nb_unique_nodes < 3) {
                continue;
            }

            // calculate centroid by averaging coordinates (uses only "real" nodes)
            size_t nb_real_nodes = 0;
            for (size_t n = 0; n < nb_cell_nodes; ++n) {
                const size_t i = size_t(cell_node_connectivity(e, n));
                if (i < firstVirtualPoint) {
                    ++nb_real_nodes;
                    centroids(e, XX) += coords(i, XX);
                    centroids(e, YY) += coords(i, YY);
                    centroids(e, ZZ) += coords(i, ZZ);
                }
            }

            if (nb_real_nodes > 1) {
                const double average_coefficient = 1. / static_cast<double>(nb_real_nodes);
                centroids(e, XX) *= average_coefficient;
                centroids(e, YY) *= average_coefficient;
                centroids(e, ZZ) *= average_coefficient;
            }
        }
    }
    return mesh.cells().field(field_name_);
}

} // namespace actions
} // namespace mesh
} // namespace atlas
