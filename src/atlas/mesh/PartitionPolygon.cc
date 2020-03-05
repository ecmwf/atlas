/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh.h"
#include "atlas/mesh/PartitionPolygon.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace mesh {

namespace {
util::Polygon::edge_set_t compute_edges( const detail::MeshImpl& mesh, idx_t halo ) {
    ATLAS_TRACE( "PartitionPolygon" );
    // extract partition boundary edges by always attempting first to`
    // remove a reversed edge from a neighbouring element, if any
    util::Polygon::edge_set_t edges;
    for ( idx_t t = 0; t < mesh.cells().nb_types(); ++t ) {
        const Elements& elements = mesh.cells().elements( t );

        const BlockConnectivity& conn = elements.node_connectivity();
        auto field_flags              = elements.view<int, 1>( elements.flags() );
        auto field_halo               = elements.view<int, 1>( elements.halo() );

        auto patch = [&field_flags]( idx_t e ) {
            using Topology = atlas::mesh::Nodes::Topology;
            return Topology::check( field_flags( e ), Topology::PATCH );
        };

        const idx_t nb_nodes = elements.nb_nodes();

        for ( idx_t j = 0; j < elements.size(); ++j ) {
            if ( patch( j ) == 0 && field_halo( j ) <= halo ) {
                for ( idx_t k = 0; k < nb_nodes; ++k ) {
                    util::Polygon::edge_t edge( conn( j, k ), conn( j, ( k + 1 ) % nb_nodes ) );
                    if ( !edges.erase( edge.reverse() ) ) {
                        edges.insert( edge );
                    }
                }
            }
        }
    }
    return edges;
}
}  // namespace

PartitionPolygon::PartitionPolygon( const detail::MeshImpl& mesh, idx_t halo ) :
    util::Polygon( compute_edges( mesh, halo ) ), mesh_( mesh ), halo_( halo ) {}

size_t PartitionPolygon::footprint() const {
    size_t size = sizeof( *this );
    size += capacity() * sizeof( idx_t );
    return size;
}

void PartitionPolygon::outputPythonScript( const eckit::PathName& filepath, const eckit::Configuration& config ) const {
    const eckit::mpi::Comm& comm = atlas::mpi::comm();
    int mpi_rank                 = int( comm.rank() );
    int mpi_size                 = int( comm.size() );

    auto xy = array::make_view<double, 2>( mesh_.nodes().xy() );

    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    for ( idx_t i : *this ) {
        xmin = std::min( xmin, xy( i, XX ) );
        xmax = std::max( xmax, xy( i, XX ) );
    }
    comm.allReduceInPlace( xmin, eckit::mpi::min() );
    comm.allReduceInPlace( xmax, eckit::mpi::max() );

    idx_t count     = mesh_.nodes().size();
    idx_t count_all = count;
    comm.allReduceInPlace( count_all, eckit::mpi::sum() );

    bool plot_nodes = config.getBool( "nodes", false );
    for ( int r = 0; r < mpi_size; ++r ) {
        if ( mpi_rank == r ) {
            std::ofstream f( filepath.asString().c_str(), mpi_rank == 0 ? std::ios::trunc : std::ios::app );

            if ( mpi_rank == 0 ) {
                f << "\n"
                     "# Configuration option to plot nodes"
                     "\n"
                     "plot_nodes = " +
                         std::string( plot_nodes ? "True" : "False" ) +
                         "\n"
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
            f << "\n"
                 "verts_"
              << r << " = [";
            for ( idx_t i : static_cast<const container_t&>( *this ) ) {
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
              << r << ".extend([Path.LINETO] * " << ( size() - 2 )
              << ")"
                 "\n"
                 "codes_"
              << r
              << ".extend([Path.CLOSEPOLY])"
                 "\n"
                 ""
                 "\n"
                 "count_"
              << r << " = " << count
              << "\n"
                 "count_all_"
              << r << " = " << count_all
              << "\n"
                 "";
            if ( plot_nodes ) {
                f << "\n"
                     "x_"
                  << r << " = [";
                for ( idx_t i = 0; i < count; ++i ) {
                    f << xy( i, XX ) << ", ";
                }
                f << "]"
                     "\n"
                     "y_"
                  << r << " = [";
                for ( idx_t i = 0; i < count; ++i ) {
                    f << xy( i, YY ) << ", ";
                }
                f << "]";
            }
            f << "\n"
                 "\n"
                 "c = cycol()"
                 "\n"
                 "ax.add_patch(patches.PathPatch(Path(verts_"
              << r << ", codes_" << r << "), facecolor=c, color=c, alpha=0.3, lw=1))";
            if ( plot_nodes ) {
                f << "\n"
                     "if plot_nodes:"
                     "\n"
                     "    ax.scatter(x_"
                  << r << ", y_" << r << ", color=c, marker='o')";
            }
            f << "\n"
                 "";
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

void PartitionPolygon::print( std::ostream& out ) const {
    out << "polygon:{"
        << "halo:" << halo_ << ",size:" << size() << ",nodes:" << static_cast<const util::Polygon&>( *this ) << "}";
}

}  // namespace mesh
}  // namespace atlas
