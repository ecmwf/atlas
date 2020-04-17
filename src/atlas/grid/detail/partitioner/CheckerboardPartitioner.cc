/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "CheckerboardPartitioner.h"

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

CheckerboardPartitioner::CheckerboardPartitioner() : Partitioner() {
    nbands_       = 0;     // to be computed later
    checkerboard_ = true;  // default
}

CheckerboardPartitioner::CheckerboardPartitioner( int N ) : Partitioner( N ) {
    nbands_       = 0;     // to be computed later
    checkerboard_ = true;  // default
}

CheckerboardPartitioner::CheckerboardPartitioner( int N, const eckit::Parametrisation & config ) : Partitioner( N ) {
    nbands_       = 0;     // to be computed later
    config.get ("nbands", nbands_);
    checkerboard_ = true;  // default
}

CheckerboardPartitioner::CheckerboardPartitioner( int N, int nbands ) : Partitioner( N ) {
    nbands_       = nbands;
    checkerboard_ = true;  // default
}

CheckerboardPartitioner::CheckerboardPartitioner( int N, int nbands, bool checkerboard ) : Partitioner( N ) {
    nbands_       = nbands;
    checkerboard_ = checkerboard;
}

CheckerboardPartitioner::Checkerboard CheckerboardPartitioner::checkerboard( const Grid& grid ) const {
    // grid dimensions
    const RegularGrid rg( grid );
    if ( !rg ) {
        throw_Exception( "Checkerboard Partitioner only works for Regular grids.", Here() );
    }

    Checkerboard cb;

    cb.nx        = rg.nx();
    cb.ny        = rg.ny();
    idx_t nparts = nb_partitions();

    if ( nbands_ > 0 ) {
        cb.nbands = nbands_;
    }
    else {
        // default number of bands
        double zz = std::sqrt( (double)( nparts * cb.ny ) / cb.nx );  // aim at +/-square regions
        cb.nbands = (idx_t)zz + 0.5;
        if ( cb.nbands < 1 ) {
            cb.nbands = 1;  // at least one band
        }
        if ( cb.nbands > nparts ) {
            cb.nbands = nparts;  // not more bands than procs
        }

        // true checkerboard means nbands must divide nparts
        if ( checkerboard_ ) {
            while ( nparts % cb.nbands != 0 ) {
                cb.nbands--;
            }
        }
    }
    if ( checkerboard_ && nparts % cb.nbands != 0 ) {
        throw_Exception( "number of bands doesn't divide number of partitions", Here() );
    }

    return cb;
}

bool compare_Y_X( const CheckerboardPartitioner::NodeInt& node1, const CheckerboardPartitioner::NodeInt& node2 ) {
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

bool compare_X_Y( const CheckerboardPartitioner::NodeInt& node1, const CheckerboardPartitioner::NodeInt& node2 ) {
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

void CheckerboardPartitioner::partition( const Checkerboard& cb, int nb_nodes, NodeInt nodes[], int part[] ) const {
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

    /*
Number of procs per band
*/
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

    /*
Number of gridpoints per band
*/
    std::vector<size_t> ngpb( nbands, 0 );
    // split latitudes?
    if ( true ) {
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

    // for (int iband=0;iband<nbands; iband++ ) std::cout << "band " << iband << "
    // : nparts = " << npartsb[iband] << "; ngpb = " << ngpb[iband] << std::endl;

    // sort nodes according to Y first, to determine bands
    std::sort( nodes, nodes + nb_nodes, compare_Y_X );

    // std::cout << __LINE__ << ",  in " << __FUNCTION__ << std::endl;

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
        // std::cout << "remainder = " << remainder << std::endl;
        for ( size_t ipart = 0; ipart < npartsb[iband]; ipart++ ) {
            ngpp[ipart] = ngpb[iband] / npartsb[iband];
            remainder -= ngpp[ipart];
            // std::cout << "remainder = " << remainder << std::endl;
        }
        // distribute remaining gridpoints over first parts
        for ( size_t ipart = 0; ipart < remainder; ipart++ ) {
            ++ngpp[ipart];
        }

        /*
std::cout << "ngpb = " << ngpb[iband] << "\n";
for (int ipart=0;ipart<npartsb[iband];ipart++) std::cout << "part " << ipart <<
"; ngpp = " << ngpp[ipart] << ". ";
std::cout << std::endl;
*/

        // set partition number for each part
        for ( size_t ipart = 0; ipart < npartsb[iband]; ipart++ ) {
            // std::cout << "looping over points " << offset << " to " <<
            // offset+ngpp[ipart] << std::endl;
            for ( size_t jj = offset; jj < offset + ngpp[ipart]; jj++ ) {
                part[nodes[jj].n] = jpart;
            }
            offset += ngpp[ipart];
            ++jpart;
        }
    }
}

void CheckerboardPartitioner::partition( const Grid& grid, int part[] ) const {
    if ( nb_partitions() == 1 )  // trivial solution, so much faster
    {
        for ( idx_t j = 0; j < grid.size(); ++j ) {
            part[j] = 0;
        }
    }
    else {
        auto cb = checkerboard( grid );

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
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::CheckerboardPartitioner>
    __CheckerBoard( "checkerboard" );
}
