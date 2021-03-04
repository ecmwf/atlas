
#include "ZonalBoardPartitioner.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <functional>
#include <iostream>
#include <vector>

#include "atlas/grid/StructuredGrid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/MicroDeg.h"

using atlas::util::microdeg;

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

ZonalBoardPartitioner::ZonalBoardPartitioner() : Partitioner() {
    nbands_       = 0;     // to be computed later
    zonalboard_ = true;  // default
}

ZonalBoardPartitioner::ZonalBoardPartitioner( int N ) : Partitioner( N ) {
    nbands_       = 0;     // to be computed later
    zonalboard_ = true;  // default
}

ZonalBoardPartitioner::ZonalBoardPartitioner( int N, const eckit::Parametrisation& config ) : Partitioner( N ) {
    config.get( "bands", nbands_ );
    config.get( "zonalboard", zonalboard_ );
}

ZonalBoardPartitioner::ZonalBoardPartitioner( int N, int nbands ) : Partitioner( N ) {
    nbands_       = nbands;
    zonalboard_ = true;  // default
}

ZonalBoardPartitioner::ZonalBoardPartitioner( int N, int nbands, bool zonalboard ) : Partitioner( N ) {
    nbands_       = nbands;
    zonalboard_ = zonalboard;
}

ZonalBoardPartitioner::Zonalboard ZonalBoardPartitioner::zonalboard( const Grid& grid ) const {
    // grid dimensions
    const RegularGrid rg( grid );
    if ( !rg ) {
        throw_Exception( "Zonalboard Partitioner only works for Regular grids.", Here() );
    }

    Zonalboard cb;

    cb.nx        = rg.nx();
    cb.ny        = rg.ny();
    idx_t nparts = nb_partitions();

    // for now try a zonal band decomposition
    cb.nbands = nparts;

    if ( zonalboard_ && nparts % cb.nbands != 0 ) {
        throw_Exception( "number of bands doesn't divide number of partitions", Here() );
    }

    return cb;
}

bool compare_Y_X( const ZonalBoardPartitioner::NodeInt& node1, const ZonalBoardPartitioner::NodeInt& node2 ) {
    // comparison of two locations; X1 < X2 if it's to the south, then to the
    // east.
    if ( node1.y < node2.y ) {
        return true;
    }
    if ( node1.y == node2.y ) {
        return ( node1.x < node2.x );
    }
    return false;
}

bool compare_X_Y( const ZonalBoardPartitioner::NodeInt& node1, const ZonalBoardPartitioner::NodeInt& node2 ) {
    // comparison of two locations; X1 < X2 if it's to the east, then to the
    // south.
    if ( node1.x < node2.x ) {
        return true;
    }
    if ( node1.x == node2.x ) {
        return ( node1.y < node2.y );
    }
    return false;
}


//  nb_nodes - total number of nodes on horizontal surface.
//  nodes[] - a C style array that holds the "global node index"?
//  part[] - a C style array that holds the partitioning number
void ZonalBoardPartitioner::partition( const Zonalboard& cb, int nb_nodes, NodeInt nodes[], int part[] ) const {
    size_t nparts = nb_partitions();
    size_t nbands = cb.nbands;
    size_t nx     = cb.nx;
    size_t ny     = cb.ny;
    size_t remainder;

    /*
Sort nodes from south to north (increasing y), and west to east (increasing x).
Now we can easily split
the points in bands. Note this may not be necessary, as it could be
already by construction in this order, but then sorting is really fast
*/
//the above comment is wrong for (at least for global) - nodes go from north to south
    /*
Number of procs per band
*/

    // ngpb - number of gridpoints per band
    std::vector<size_t> ngpb( nbands, 0 );

    remainder = ny;
    // syslice - standard y slice for all partitions except South Pole that holds the remainder.
    std::size_t syslice = ny / nparts;
    for ( std::size_t iband = 0; iband < nbands -1; iband++ ) {
      ngpb[iband] = nx * syslice;
      remainder -= ngpb[iband] / nx;
    }
    ngpb[nbands -1] = remainder * nx;

    for (int iband=0;iband<nbands; iband++ ) std::cout << "band " << iband <<
      "; ngpb = " << ngpb[iband] << std::endl;

    // sort nodes according to Y first, to determine bands
    std::sort( nodes, nodes + nb_nodes, compare_Y_X );

    // std::cout << __LINE__ << ",  in " << __FUNCTION__ << std::endl;

    // for each band, select gridpoints belonging to that band, and sort them
    // according to X first
    size_t offset = 0;
    for ( size_t iband = 0; iband < nbands; iband++ ) {

      // sort according to X first
      std::sort( nodes + offset, nodes + offset + ngpb[iband], compare_X_Y );

      int ngpp(ngpb[iband]);
      for ( size_t jj = offset; jj < offset + ngpb[iband]; jj++ ) {
          part[nodes[jj].n] = iband;
      }
      offset += ngpp;
    }
}

void ZonalBoardPartitioner::partition( const Grid& grid, int part[] ) const {
    if ( nb_partitions() == 1 )  // trivial solution, so much faster
    {
        for ( idx_t j = 0; j < grid.size(); ++j ) {
            part[j] = 0;
        }
    }
    else {
        auto cb = zonalboard( grid );

        std::vector<NodeInt> nodes( grid.size() );
        int n( 0 );

        for ( idx_t iy = 0; iy < cb.ny; ++iy ) {
            for ( idx_t ix = 0; ix < cb.nx; ++ix ) {
                nodes[n].x = static_cast<int>( ix );
                nodes[n].y = static_cast<int>( iy );
                nodes[n].n = static_cast<int>( n );
                ++n;
            }
        }

        partition( cb, grid.size(), nodes.data(), part );
    }
}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas

namespace {
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::ZonalBoardPartitioner>
    __ZonalBoard( "zonalboard" );
}
