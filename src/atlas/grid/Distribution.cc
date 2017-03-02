/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include "atlas/grid/Grid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/detail/partitioners/Partitioner.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace grid {

Distribution::impl_t::impl_t( const Grid& grid ) :
    nb_partitions_(1),
    part_(grid.npts(),0),
    nb_pts_(nb_partitions_,grid.npts()),
    max_pts_(grid.npts()),
    min_pts_(grid.npts()) {
}

Distribution::impl_t::impl_t( const Partitioner& partitioner ) {
    part_.resize(partitioner.grid().npts());
    partitioner.partition(part_.data());
    nb_partitions_ = partitioner.nb_partitions();
    nb_pts_.resize(nb_partitions_,0);
    for(size_t j = 0; j < part_.size(); ++j)
        ++nb_pts_[part_[j]];
    max_pts_ = *std::max_element(nb_pts_.begin(),nb_pts_.end());
    min_pts_ = *std::min_element(nb_pts_.begin(),nb_pts_.end());
}

Distribution::impl_t::impl_t( size_t npts, int part[], int part0 ) {
    part_.assign(part,part+npts);
    std::set<int> partset(part_.begin(),part_.end());
    nb_partitions_ = partset.size();
    nb_pts_.resize(nb_partitions_,0);
    for(size_t j = 0; j < part_.size(); ++j) {
        part_[j] -= part0;
        ++nb_pts_[part_[j]];
    }
    max_pts_ = *std::max_element(nb_pts_.begin(),nb_pts_.end());
    min_pts_ = *std::min_element(nb_pts_.begin(),nb_pts_.end());
}



Distribution::Distribution():
    impl_( nullptr ) {
}

Distribution::Distribution( const impl_t* impl ):
    impl_( impl ) {
}

Distribution::Distribution( const Distribution& other ):
    impl_( other.impl_ ) {
}

Distribution::Distribution( const Grid& grid ):
    impl_( new impl_t(grid) ) {
}

Distribution::Distribution(const Partitioner& partitioner):
    impl_( new impl_t(partitioner) ) {
}

Distribution::Distribution(size_t npts, int part[], int part0):
    impl_( new impl_t(npts,part,part0) ) {
}



Distribution::impl_t* atlas__GridDistribution__new(int npts, int part[], int part0) {
    Distribution::impl_t* dist = new Distribution::impl_t(npts,part,part0);
    Log::info() << Here() << "owners = " << dist->owners() << std::endl;

    return dist;
//    return new GridDistribution::impl_t(npts,part,part0);
}

void atlas__GridDistribution__delete(Distribution::impl_t* This) {
    Log::info() << Here() << "owners = " << This->owners() << std::endl;
    delete This;
}

} // namespace grid
} // namespace atlas
