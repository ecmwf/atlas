/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <sstream>

#include "atlas/array.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/trans/ifs/TransIFS.h"

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

TransPartitioner::TransPartitioner() : Partitioner() {
    EqualRegionsPartitioner eqreg( nb_partitions() );
    nbands_ = eqreg.nb_bands();
    nregions_.resize( nbands_ );
    for ( size_t b = 0; b < nbands_; ++b ) {
        nregions_[b] = eqreg.nb_regions( b );
    }
}

TransPartitioner::TransPartitioner( const idx_t N, const eckit::Parametrisation& ) : Partitioner( N ) {
    EqualRegionsPartitioner eqreg( nb_partitions() );
    nbands_ = eqreg.nb_bands();
    nregions_.resize( nbands_ );
    for ( size_t b = 0; b < nbands_; ++b ) {
        nregions_[b] = eqreg.nb_regions( b );
    }
}

TransPartitioner::~TransPartitioner() = default;

void TransPartitioner::partition( const Grid& grid, int part[] ) const {
    ATLAS_TRACE( "TransPartitioner::partition" );

    StructuredGrid g( grid );
    if ( not g ) {
        throw_Exception( "Grid is not a grid::Structured type. Cannot partition using IFS trans", Here() );
    }

    trans::TransIFS t( grid );
    if ( nb_partitions() != idx_t( t.nproc() ) ) {
        std::stringstream msg;
        msg << "Requested to partition grid with TransPartitioner in " << nb_partitions()
            << " partitions, but "
               "the internal Trans library could only be set-up for "
            << t.nproc()
            << " partitions "
               "(equal to number of MPI tasks in communicator).";
        throw_Exception( msg.str(), Here() );
    }

    int nlonmax = g.nxmax();

    array::LocalView<int, 1> nloen       = t.nloen();
    array::LocalView<int, 1> n_regions   = t.n_regions();
    array::LocalView<int, 1> nfrstlat    = t.nfrstlat();
    array::LocalView<int, 1> nlstlat     = t.nlstlat();
    array::LocalView<int, 1> nptrfrstlat = t.nptrfrstlat();
    array::LocalView<int, 2> nsta        = t.nsta();
    array::LocalView<int, 2> nonl        = t.nonl();

    int i( 0 );
    int maxind( 0 );
    std::vector<int> iglobal( t.ndgl() * nlonmax, -1 );

    for ( int jgl = 0; jgl < t.ndgl(); ++jgl ) {
        for ( int jl = 0; jl < nloen( jgl ); ++jl ) {
            ++i;
            iglobal[jgl * nlonmax + jl] = i;
            maxind                      = std::max( maxind, jgl * nlonmax + jl );
        }
    }
    int iproc( 0 );
    for ( int ja = 0; ja < t.n_regions_NS(); ++ja ) {
        for ( int jb = 0; jb < n_regions( ja ); ++jb ) {
            for ( int jgl = nfrstlat( ja ) - 1; jgl < nlstlat( ja ); ++jgl ) {
                int igl = nptrfrstlat( ja ) + jgl - nfrstlat( ja );
                for ( int jl = nsta( jb, igl ) - 1; jl < nsta( jb, igl ) + nonl( jb, igl ) - 1; ++jl ) {
                    idx_t ind = iglobal[jgl * nlonmax + jl] - 1;
                    if ( ind >= grid.size() ) {
                        throw_OutOfRange( "part", ind, grid.size(), Here() );
                    }
                    part[ind] = iproc;
                }
            }
            ++iproc;
        }
    }
}

int TransPartitioner::nb_bands() const {
    return nbands_;
}

int TransPartitioner::nb_regions( int b ) const {
    return nregions_[b];
}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas

namespace {
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::TransPartitioner> __Trans(
    "trans" );
atlas::grid::detail::partitioner::PartitionerBuilder<atlas::grid::detail::partitioner::TransPartitioner> __TransIFS(
    "ifs" );
}  // namespace
