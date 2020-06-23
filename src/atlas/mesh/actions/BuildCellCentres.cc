/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <limits>

#include "eckit/types/FloatCompare.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

BuildCellCentres::BuildCellCentres( const std::string& field_name, bool force_recompute ) :
    field_name_( field_name ), force_recompute_( force_recompute ), flatten_virtual_elements_( true ) {}

BuildCellCentres::BuildCellCentres( eckit::Configuration& config ) :
    field_name_( config.getString( "name", "centre" ) ),
    force_recompute_( config.getBool( "force_recompute", false ) ),
    flatten_virtual_elements_( config.getBool( "flatten_virtual_elements", true ) ) {}

Field& BuildCellCentres::operator()( Mesh& mesh ) const {
    bool recompute = force_recompute_;
    if ( !mesh.cells().has_field( field_name_ ) ) {
        mesh.cells().add(
            Field( field_name_, array::make_datatype<double>(), array::make_shape( mesh.cells().size(), 3 ) ) );
        recompute = true;
    }
    if ( recompute ) {
        ATLAS_TRACE( "BuildCellCentres" );
        mesh::Nodes& nodes                 = mesh.nodes();
        array::ArrayView<double, 2> coords = array::make_view<double, 2>( nodes.field( "xyz" ) );

        idx_t firstVirtualPoint = std::numeric_limits<idx_t>::max();
        if ( nodes.metadata().has( "NbRealPts" ) ) {
            firstVirtualPoint = nodes.metadata().get<idx_t>( "NbRealPts" );
        }

        idx_t nb_cells = mesh.cells().size();
        auto centroids = array::make_view<double, 2>( mesh.cells().field( field_name_ ) );
        const mesh::HybridElements::Connectivity& cell_node_connectivity = mesh.cells().node_connectivity();

        for ( idx_t e = 0; e < nb_cells; ++e ) {
            centroids( e, XX ) = 0.;
            centroids( e, YY ) = 0.;
            centroids( e, ZZ ) = 0.;

            const idx_t nb_cell_nodes = cell_node_connectivity.cols( e );

            // check for degenerate elements (less than three unique nodes)
            // NOTE: this is not a proper check but it is very robust
            eckit::types::CompareApproximatelyEqual<double> approx( 1.e-9 );

            idx_t nb_equal_nodes = 0;
            for ( idx_t ni = 0; ni < nb_cell_nodes - 1; ++ni ) {
                idx_t i = cell_node_connectivity( e, ni );
                Point3 Pi( coords( i, XX ), coords( i, YY ), coords( i, ZZ ) );
                for ( idx_t nj = ni + 1; nj < nb_cell_nodes; ++nj ) {
                    idx_t j = cell_node_connectivity( e, nj );
                    Point3 Pj( coords( j, XX ), coords( j, YY ), coords( j, ZZ ) );
                    if ( approx( Pi[XX], Pj[XX] ) && approx( Pi[YY], Pj[YY] ) && approx( Pi[ZZ], Pj[ZZ] ) ) {
                        ++nb_equal_nodes;
                    }
                }
            }

            idx_t nb_unique_nodes = idx_t( nb_cell_nodes ) - nb_equal_nodes;
            if ( nb_unique_nodes < 3 ) {
                continue;
            }

            if ( flatten_virtual_elements_ ) {
                // calculate centroid by averaging coordinates (uses only "real" nodes)
                idx_t nb_real_nodes = 0;
                for ( idx_t n = 0; n < nb_cell_nodes; ++n ) {
                    const idx_t i = cell_node_connectivity( e, n );
                    if ( i < firstVirtualPoint ) {
                        ++nb_real_nodes;
                        centroids( e, XX ) += coords( i, XX );
                        centroids( e, YY ) += coords( i, YY );
                        centroids( e, ZZ ) += coords( i, ZZ );
                    }
                }

                if ( nb_real_nodes > 1 ) {
                    const double average_coefficient = 1. / static_cast<double>( nb_real_nodes );
                    centroids( e, XX ) *= average_coefficient;
                    centroids( e, YY ) *= average_coefficient;
                    centroids( e, ZZ ) *= average_coefficient;
                }
            }
            else {
                const double average_coefficient = 1. / static_cast<double>( nb_cell_nodes );
                for ( idx_t n = 0; n < nb_cell_nodes; ++n ) {
                    const idx_t i = cell_node_connectivity( e, n );
                    for ( idx_t d = 0; d < 3; ++d ) {
                        centroids( e, d ) += coords( i, d ) * average_coefficient;
                    }
                }
            }
        }
    }
    return mesh.cells().field( field_name_ );
}

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
