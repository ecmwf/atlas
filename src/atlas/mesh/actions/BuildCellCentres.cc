/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/field/Field.h"
#include "atlas/internals/Parameters.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

void BuildCellCentres::operator()( Mesh& mesh ) const
{
    if( mesh.cells().has_field("centre") ) {

        mesh::Nodes& nodes     = mesh.nodes();
        array::ArrayView<double,2> coords  ( nodes.field("xyz") );

        size_t nb_cells = mesh.cells().size();
        array::ArrayView<double,2> centroids ( mesh.cells().add( field::Field::create<double>("centre", array::make_shape(nb_cells,3))) );
        const mesh::HybridElements::Connectivity& cell_node_connectivity = mesh.cells().node_connectivity();

        for (size_t e=0; e<nb_cells; ++e)
        {
            centroids(e,internals::XX) = 0.;
            centroids(e,internals::YY) = 0.;
            centroids(e,internals::ZZ) = 0.;
            const size_t nb_nodes_per_elem = cell_node_connectivity.cols(e);
            const double average_coefficient = 1./static_cast<double>(nb_nodes_per_elem);
            for (size_t n=0; n<nb_nodes_per_elem; ++n)
            {
                centroids(e,internals::XX) += coords( cell_node_connectivity(e,n), internals::XX );
                centroids(e,internals::YY) += coords( cell_node_connectivity(e,n), internals::YY );
                centroids(e,internals::ZZ) += coords( cell_node_connectivity(e,n), internals::ZZ );
            }
            centroids(e,internals::XX) *= average_coefficient;
            centroids(e,internals::YY) *= average_coefficient;
            centroids(e,internals::ZZ) *= average_coefficient;
        }
    }
}

} // namespace actions
} // namespace mesh
} // namespace atlas
