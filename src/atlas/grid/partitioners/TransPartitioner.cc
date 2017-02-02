/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/exception/Exceptions.h"
#include "atlas/grid/Structured.h"
#include "atlas/grid/partitioners/TransPartitioner.h"
#include "atlas/trans/Trans.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace grid {
namespace partitioners {

TransPartitioner::TransPartitioner(const grid::Grid& grid,
                                   const trans::Trans& trans) :
    Partitioner(grid, trans.nproc()),
    t_( const_cast<trans::Trans*>(&trans)), owned_(false) {
    ASSERT( t_ != NULL );
}

TransPartitioner::TransPartitioner(const grid::Grid& grid) :
    Partitioner(grid),
    t_( new trans::Trans(grid,0) ), owned_(true) {
    ASSERT( t_ != NULL );
    ASSERT( size_t(t_->nproc()) == nb_partitions() );
}

TransPartitioner::TransPartitioner(const grid::Grid& grid,
                                   const size_t nb_partitions) :
    Partitioner(grid, nb_partitions),
    t_( new trans::Trans(grid,0) ), owned_(true) {
    ASSERT( t_ != NULL );
    if( nb_partitions != size_t(t_->nproc()) ) {
        std::stringstream msg;
        msg << "Requested to partition grid with TransPartitioner in "<<nb_partitions<<" partitions, but "
            "the internal Trans library could only be set-up for "<<t_->nproc()<< " partitions "
            "(equal to number of MPI tasks in communicator).";
        throw eckit::BadParameter(msg.str(),Here());
    }
}

TransPartitioner::~TransPartitioner() {
    if( owned_ )
        delete t_;
}

void TransPartitioner::partition(int part[]) const {
    if( dynamic_cast<const grid::Structured*>(&grid()) == NULL )
        throw eckit::BadCast("Grid is not a grid::Structured type. Cannot partition using IFS trans",Here());

    int nlonmax = dynamic_cast<const grid::Structured*>(&grid())->nlonmax();

    array::LocalView<int,1> nloen       = t_->nloen();
    array::LocalView<int,1> n_regions   = t_->n_regions();
    array::LocalView<int,1> nfrstlat    = t_->nfrstlat();
    array::LocalView<int,1> nlstlat     = t_->nlstlat();
    array::LocalView<int,1> nptrfrstlat = t_->nptrfrstlat();
    array::LocalView<int,2> nsta        = t_->nsta();
    array::LocalView<int,2> nonl        = t_->nonl();


    int i(0);
    int maxind(0);
    std::vector<int> iglobal(t_->ndgl()*nlonmax,-1);

    for( int jgl=0; jgl<t_->ndgl(); ++jgl ) {
        for( int jl=0; jl<nloen(jgl); ++jl ) {
            ++i;
            iglobal[jgl*nlonmax+jl] = i;
            maxind = std::max(maxind,jgl*nlonmax+jl);
        }
    }
    int iproc(0);
    for( int ja=0; ja<t_->n_regions_NS(); ++ja ) {
        for( int jb=0; jb<n_regions(ja); ++jb ) {
            for( int jgl=nfrstlat(ja)-1; jgl<nlstlat(ja); ++jgl ) {
                int igl = nptrfrstlat(ja) + jgl - nfrstlat(ja);
                for( int jl=nsta(jb,igl)-1; jl<nsta(jb,igl)+nonl(jb,igl)-1; ++jl ) {
                    size_t ind = iglobal[jgl*nlonmax+jl] - 1;
                    if( ind >= grid().npts() ) throw eckit::OutOfRange(ind,grid().npts(),Here());
                    part[ind] = iproc;
                }
            }
            ++iproc;
        }
    }
}

int TransPartitioner::nb_bands() const {
    ASSERT( t_!= NULL );
    return t_->n_regions_NS();
}

int TransPartitioner::nb_regions(int b) const {
    ASSERT( t_!= NULL );
    return t_->n_regions()(b);
}

} // namespace partitioners
} // namespace grid
} // namespace atlas

namespace {
atlas::grid::partitioners::PartitionerBuilder<
atlas::grid::partitioners::TransPartitioner> __Trans("Trans");
}

