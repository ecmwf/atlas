/*
 * (C) Copyright 1996-2017 ECMWF.
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
#include "atlas/util/CoordinateEnums.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

void BuildCellCentres::operator()( Mesh& mesh ) const
{
    if( !mesh.cells().has_field("centre") ) {

        mesh::Nodes& nodes     = mesh.nodes();
        array::ArrayView<double,2> coords = array::make_view<double,2>( nodes.field("xyz") );

        size_t nb_cells = mesh.cells().size();
        array::ArrayView<double,2> centroids = array::make_view<double,2>( 
          mesh.cells().add(
           field::Field("centre", array::make_datatype<double>(), array::make_shape(nb_cells,3))) );
        const mesh::HybridElements::Connectivity& cell_node_connectivity = mesh.cells().node_connectivity();

        for (size_t e=0; e<nb_cells; ++e)
        {
            centroids(e,XX) = 0.;
            centroids(e,YY) = 0.;
            centroids(e,ZZ) = 0.;
            const size_t nb_nodes_per_elem = cell_node_connectivity.cols(e);
            const double average_coefficient = 1./static_cast<double>(nb_nodes_per_elem);
            for (size_t n=0; n<nb_nodes_per_elem; ++n)
            {
                centroids(e,XX) += coords( cell_node_connectivity(e,n), XX );
                centroids(e,YY) += coords( cell_node_connectivity(e,n), YY );
                centroids(e,ZZ) += coords( cell_node_connectivity(e,n), ZZ );
            }
            centroids(e,XX) *= average_coefficient;
            centroids(e,YY) *= average_coefficient;
            centroids(e,ZZ) *= average_coefficient;
        }
    }
}

} // namespace actions
} // namespace mesh
} // namespace atlas
