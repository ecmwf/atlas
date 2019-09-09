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

#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/ReorderHilbert.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

void outputHilbertCurve( const eckit::PathName& filepath, const array::ArrayView<double, 2>& xy ) {
    const eckit::mpi::Comm& comm = atlas::mpi::comm();
    int mpi_rank                 = int( comm.rank() );
    int mpi_size                 = int( comm.size() );

    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    for ( idx_t i = 0; i < xy.shape( 0 ); ++i ) {
        xmin = std::min( xmin, xy( i, XX ) );
        xmax = std::max( xmax, xy( i, XX ) );
    }
    comm.allReduceInPlace( xmin, eckit::mpi::min() );
    comm.allReduceInPlace( xmax, eckit::mpi::max() );

    idx_t count     = xy.shape( 0 );
    idx_t count_all = count;
    comm.allReduceInPlace( count_all, eckit::mpi::sum() );

    for ( int r = 0; r < mpi_size; ++r ) {
        if ( mpi_rank == r ) {
            std::ofstream f( filepath.asString().c_str(), mpi_rank == 0 ? std::ios::trunc : std::ios::app );

            if ( mpi_rank == 0 ) {
                f << "\n"
                     "\n"
                     "import matplotlib.pyplot as plt"
                     "\n"
                     "from matplotlib.path import Path"
                     "\n"
                     "import matplotlib.patches as patches"
                     "\n"
                     ""
                     "\n"
                     "from itertools import cycle"
                     "\n"
                     "import matplotlib.cm as cm"
                     "\n"
                     "import numpy as np"
                     "\n"
                     "cycol = cycle([cm.Paired(i) for i in "
                     "np.linspace(0,1,12,endpoint=True)]).next"
                     "\n"
                     ""
                     "\n"
                     "fig = plt.figure()"
                     "\n"
                     "ax = fig.add_subplot(111,aspect='equal')"
                     "\n"
                     "";
            }

            if ( mpi_rank == r ) {  // replace "r" with rank you wish to plot only
                f << "\n"
                     "verts_"
                  << r << " = [";
                for ( idx_t n = 0; n < xy.shape( 0 ); ++n ) {
                    idx_t i = n;  //reorder.hilbert_reordering_[n].second;
                    f << "\n  (" << xy( i, XX ) << ", " << xy( i, YY ) << "), ";
                }
                f << "\n]"
                     "\n"
                     ""
                     "\n"
                     "codes_"
                  << r
                  << " = [Path.MOVETO]"
                     "\n"
                     "codes_"
                  << r << ".extend([Path.LINETO] * " << ( xy.shape( 0 ) - 2 )
                  << ")"
                     "\n"
                     "codes_"
                  << r
                  << ".extend([Path.LINETO])"
                     "\n"
                     ""
                     "\n"
                     "count_"
                  << r << " = " << count
                  << "\n"
                     "count_all_"
                  << r << " = " << count_all
                  << "\n"
                     "c = cycol()"
                     "\n"
                     "xs_"
                  << r << ", ys_" << r << " = zip(*verts_" << r
                  << ")"
                     "\n"
                     "ax.plot(xs_"
                  << r << ",ys_" << r
                  << ", '-', lw=1, color=c )"
                     "\n"
                     //                 "for i in range( len(verts_0) ):"
                     //                 "\n"
                     //                 "  plt.text( xs_0[i], ys_0[i], str(i) )"
                     //                 "\n"
                     "";
            }
            if ( mpi_rank == mpi_size - 1 ) {
                f << "\n"
                     "ax.set_xlim( "
                  << xmin << "-5, " << xmax
                  << "+5)"
                     "\n"
                     "ax.set_ylim(-90-5,  90+5)"
                     "\n"
                     "ax.set_xticks([0,45,90,135,180,225,270,315,360])"
                     "\n"
                     "ax.set_yticks([-90,-45,0,45,90])"
                     "\n"
                     "plt.grid()"
                     "\n"
                     "plt.show()";
            }
        }
        comm.barrier();
    }
}

Field create_cell_centres( Mesh& mesh ) {
    auto cell_centres =
        Field( "cell_centres", array::make_datatype<double>(), array::make_shape( mesh.cells().size(), 2 ) );
    auto nodes_xy = array::make_view<double, 2>( mesh.nodes().xy() );
    for ( idx_t t = 0; t < mesh.cells().nb_types(); ++t ) {
        auto& cells = mesh.cells().elements( t );
        auto xy     = cells.view<double, 2>( cell_centres );

        // Compute cell-centres
        {
            const auto& node_connectivity = cells.node_connectivity();
            const idx_t nb_nodes          = cells.nb_nodes();
            const double nb_nodes_double  = nb_nodes;
            for ( idx_t e = 0; e < cells.size(); ++e ) {
                double x{0};
                double y{0};
                for ( idx_t c = 0; c < nb_nodes; ++c ) {
                    idx_t n = node_connectivity( e, c );
                    x += nodes_xy( n, XX );
                    y += nodes_xy( n, YY );
                }
                xy( e, XX ) = x / nb_nodes_double;
                xy( e, YY ) = y / nb_nodes_double;
            }
        }
    }
    return cell_centres;
}

Field create_edge_centres( Mesh& mesh ) {
    auto edge_centres =
        Field( "edge_centres", array::make_datatype<double>(), array::make_shape( mesh.edges().size(), 2 ) );
    auto nodes_xy = array::make_view<double, 2>( mesh.nodes().xy() );
    auto& edges   = mesh.edges();
    auto xy       = array::make_view<double, 2>( edge_centres );

    // Compute edge-centres
    {
        const auto& node_connectivity = edges.node_connectivity();
        for ( idx_t e = 0; e < edges.size(); ++e ) {
            double x{0};
            double y{0};
            for ( idx_t c = 0; c < 2; ++c ) {
                idx_t n = node_connectivity( e, c );
                x += nodes_xy( n, XX );
                y += nodes_xy( n, YY );
            }
            xy( e, XX ) = 0.5 * x;
            xy( e, YY ) = 0.5 * y;
        }
    }
    return edge_centres;
}

CASE( "test_hilbert_reordering" ) {
    auto meshgenerator = StructuredMeshGenerator( util::Config( "patch_pole", false )( "triangulate", true ) );
    Mesh mesh          = meshgenerator( Grid( "O32" ) );
    auto reorder       = mesh::actions::ReorderHilbert{util::Config( "recursion", 30 )};
    reorder( mesh );

    // functionspace::NodeColumns( mesh, option::halo( 2 ) );

    auto xy = array::make_view<double, 2>( mesh.nodes().xy() );
    outputHilbertCurve( "hilbert_nodes.py", xy );
    output::Gmsh gmsh{"hilbert.msh", util::Config( "coordinates", "xy" )};
    gmsh.write( mesh );

    Field cell_centres   = create_cell_centres( mesh );
    auto cell_centres_xy = array::make_view<double, 2>( cell_centres );
    outputHilbertCurve( "hilbert_elements.py", cell_centres_xy );

    mesh::actions::build_edges( mesh );
    Field edge_centres   = create_edge_centres( mesh );
    auto edge_centres_xy = array::make_view<double, 2>( edge_centres );
    outputHilbertCurve( "hilbert_edges.py", edge_centres_xy );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
