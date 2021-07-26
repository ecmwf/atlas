/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "CubedSpherePartitioner.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <functional>
#include <iostream>
#include <vector>


#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/MicroDeg.h"

using atlas::util::microdeg;

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

namespace  {

bool isNearInt( double value )
{
    const double epsilon = 1e-12;
    const double diff = value - floor( value );
    return ( diff <= epsilon || diff >= (1.0 - epsilon) );
}

}

CubedSpherePartitioner::CubedSpherePartitioner() : Partitioner() {}

CubedSpherePartitioner::CubedSpherePartitioner( int N ) : Partitioner( N ) {}

CubedSpherePartitioner::CubedSpherePartitioner( int N, const eckit::Parametrisation& config ) : Partitioner( N ) {
    config.get( "starting rank on tile", globalProcStartPE_);
    config.get( "final rank on tile", globalProcEndPE_);
    config.get( "nprocy", nprocy_);
    config.get( "nprocx", nprocx_);
    config.get( "nprocy", nprocy_);
}

CubedSpherePartitioner::CubedSpherePartitioner( const int N, const std::vector<int> & globalProcStartPE,
                                                const std::vector<int> & globalProcEndPE,
                                                const std::vector<int> & nprocx,
                                                const std::vector<int> & nprocy )
    : Partitioner( N ), globalProcStartPE_{globalProcStartPE},  globalProcEndPE_{globalProcEndPE},
      nprocx_{nprocx},  nprocy_{nprocy} {}



CubedSpherePartitioner::CubedSphere CubedSpherePartitioner::cubedsphere( const Grid& grid ) const {
    // grid dimensions
    const CubedSphereGrid cg( grid );
    if ( !cg ) {
        throw_Exception( "CubedSphere Partitioner only works for Regular grids.", Here() );
    }

    CubedSphere cb;

    for (idx_t t = 0; t < 6; ++t) {
        cb.nx[t] = cg.N();
        cb.ny[t] = cg.N();
    }

    atlas::idx_t nparts = nb_partitions();

    if (regular_) {
      // share PEs around tiles
      // minRanksPerTile

      idx_t ranksPerTile = nparts/6;

      idx_t reminder = nparts - 6 * ranksPerTile;

      for (idx_t t = 0; t < 6; ++t) {
         cb.nproc[t] = ranksPerTile;
      }

      // round-robin;
      idx_t t{0};
      while (reminder > 0) {
         if (t == 6) t=0;
         cb.nproc[t] += 1;
         t += 1;
         reminder -= 1;
      }

      // now need to specify nprocx and nprocy.
      // nproc is 0 for the tile we use default nprocx and nprocy = 1

      // if we can square-root nproc and get an integer
      // we use that for nprocx and nprocy
      // otherwise we split just in nprocx and keep nprocy =1;

       for (idx_t t = 0; t < 6; ++t) {
           if (cb.nproc[t] > 0) {
              double sq = std::sqrt(static_cast<double>(cb.nproc[t]));
           }
       }
    }
    return cb;
}

bool compare_Y_X( const CubedSpherePartitioner::NodeInt& node1, const CubedSpherePartitioner::NodeInt& node2 ) {
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

bool compare_X_Y( const CubedSpherePartitioner::NodeInt& node1, const CubedSpherePartitioner::NodeInt& node2 ) {
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

void CubedSpherePartitioner::partition( const CubedSphere& cb, const int nb_nodes, NodeInt nodes[], int part[] ) const {
    size_t nparts = nb_partitions();
    long remainder;

    /* Sort cell centers (increasing tiles), south to north (increasing y) and west to east (increasing x)

    /* No of procs per tile

    /*


    /*
Sort nodes from south to north (increasing y), and west to east (increasing x).
Now we can easily split
the points in bands. Note this may not be necessary, as it could be
already by construction in this order, but then sorting is really fast
*/

    /*
Number of procs per band

    std::vector<size_t> npartsb( nbands, 0 );  // number of procs per band
    remainder = nparts;
    for ( size_t iband = 0; iband < nbands; iband++ ) {
        npartsb[iband] = nparts / nbands;
        remainder -= npartsb[iband];
    }
    // distribute remaining procs over first bands
    for ( size_t iband = 0; iband < remainder; iband++ ) {
        ++npartsb[iband];
    }

    bool split_lons = not regular_;
    bool split_lats = not regular_;

*/

    /*
Number of gridpoints per band
*/


    /*
    std::vector<size_t> ngpb( nbands, 0 );
    // split latitudes?
    if ( split_lats ) {
        remainder = nb_nodes;
        for ( size_t iband = 0; iband < nbands; iband++ ) {
            ngpb[iband] = ( nb_nodes * npartsb[iband] ) / nparts;
            remainder -= ngpb[iband];
        }
        // distribute remaining gridpoints over first bands
        for ( size_t iband = 0; iband < remainder; iband++ ) {
            ++ngpb[iband];
        }
    }
    else {
        remainder = ny;
        for ( size_t iband = 0; iband < nbands; iband++ ) {
            ngpb[iband] = nx * ( ( ny * npartsb[iband] ) / nparts );
            remainder -= ngpb[iband] / nx;
        }
        // distribute remaining rows over first bands
        for ( size_t iband = 0; iband < remainder; iband++ ) {
            ngpb[iband] += nx;
        }
    }

    // sort nodes according to Y first, to determine bands
    std::sort( nodes, nodes + nb_nodes, compare_Y_X );

    // for each band, select gridpoints belonging to that band, and sort them
    // according to X first
    size_t offset = 0;
    int jpart     = 0;
    for ( size_t iband = 0; iband < nbands; iband++ ) {
        // sort according to X first
        std::sort( nodes + offset, nodes + offset + ngpb[iband], compare_X_Y );

        // number of gridpoints per task
        std::vector<int> ngpp( npartsb[iband], 0 );
        remainder = ngpb[iband];

        int part_ny = ngpb[iband] / cb.nx;
        int part_nx = ngpb[iband] / npartsb[iband] / part_ny;

        for ( size_t ipart = 0; ipart < npartsb[iband]; ipart++ ) {
            if ( split_lons ) {
                ngpp[ipart] = ngpb[iband] / npartsb[iband];
            }
            else {
                ngpp[ipart] = part_nx * part_ny;
            }
            remainder -= ngpp[ipart];
        }
        if ( split_lons ) {
            // distribute remaining gridpoints over first parts
            for ( size_t ipart = 0; ipart < remainder; ipart++ ) {
                ++ngpp[ipart];
            }
        }
        else {
            size_t ipart = 0;
            while ( remainder > part_ny ) {
                ngpp[ipart++] += part_ny;
                remainder -= part_ny;
            }
            ngpp[npartsb[iband] - 1] += remainder;
        }

        // set partition number for each part
        for ( size_t ipart = 0; ipart < npartsb[iband]; ipart++ ) {
            for ( size_t jj = offset; jj < offset + ngpp[ipart]; jj++ ) {
                part[nodes[jj].n] = jpart;
            }
            offset += ngpp[ipart];
            ++jpart;
        }
    }
    */
}

void CubedSpherePartitioner::partition( const Grid& grid, int part[] ) const {
    if ( nb_partitions() == 1 )  // trivial solution, so much faster
    {
        for ( idx_t j = 0; j < grid.size(); ++j ) {
            part[j] = 0;
        }
    }
    else {
        auto cb = cubedsphere( grid );

        std::vector<NodeInt> nodes( grid.size() );
        int n( 0 );

        for (idx_t it = 0; it < 6; ++it) {
            for ( idx_t iy = 0; iy < cb.ny[it]; ++iy ) {
                for ( idx_t ix = 0; ix < cb.nx[it]; ++ix ) {
                    nodes[n].t = static_cast<int>( it );
                    nodes[n].x = static_cast<int>( ix );
                    nodes[n].y = static_cast<int>( iy );
                    nodes[n].n = static_cast<int>( n );
                    ++n;
                }
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
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::CubedSpherePartitioner>
    __CubedSphere( "CubedSphere" );
}
