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

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

#include "atlas/mesh/actions/ReorderHilbert.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

void outputPythonScript( const eckit::PathName& filepath, const eckit::Configuration& config,
                         const array::ArrayView<double, 2>& xy, const mesh::actions::ReorderHilbert& reorder ) {
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
                     "# Configuration option to plot nodes"
                     "\n"
                     "plot_nodes = " +
                         std::string( ( config.getBool( "nodes", false ) ? "True" : "False" ) ) +
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
                 ""
                 //"\n" "x_" << r << " = ["; for (idx_t i=0; i<count; ++i) { f <<
                 // xy(i, XX) << ", "; } f << "]"
                 //"\n" "y_" << r << " = ["; for (idx_t i=0; i<count; ++i) { f <<
                 // xy(i, YY) << ", "; } f << "]"
                 "\n"
                 "\n"
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


CASE( "test_hilbert_reordering" ) {
    Mesh mesh    = StructuredMeshGenerator().generate( Grid( "O32" ) );
    auto reorder = mesh::actions::ReorderHilbert{mesh, util::Config( "recursion", 8 )};
    auto xy      = array::make_view<double, 2>( mesh.nodes().xy() );

    reorder();

#if 1
    gidx_t previous_hilbert_idx{-1};
    for ( idx_t n = 0; n < xy.shape( 0 ); ++n ) {
        auto hilbert_idx = reorder.hilbert_reordering_[n].first;
        idx_t ip         = reorder.hilbert_reordering_[n].second;
        PointXY p{xy( ip, 0 ), xy( ip, 1 )};
        //Log::info() << ip << '\t' << p << "\t" << hilbert_idx << std::endl;
        if ( hilbert_idx == previous_hilbert_idx ) {
            //throw_Exception("Duplicate hilbert index detected in ReorderHilbert", Here() );
        }
        previous_hilbert_idx = hilbert_idx;
    }

    outputPythonScript( "hilbert.py", util::NoConfig(), xy, reorder );
#endif

    output::Gmsh gmsh{"hilbert.msh"};
    gmsh.write( mesh );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
