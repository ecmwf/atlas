/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>

#include "Distribution.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/distribution/DistributionImpl.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace grid {

Distribution::Distribution() : impl_( nullptr ) {}

Distribution::Distribution( const Implementation* impl ) : impl_( impl ) {}

Distribution::Distribution( const Distribution& other ) : impl_( other.impl_ ) {}

Distribution::Distribution( const Grid& grid ) : impl_( new Implementation( grid ) ) {}

Distribution::Distribution( const Grid& grid, const Partitioner& partitioner ) :
    impl_( new Implementation( grid, partitioner ) ) {}

Distribution::Distribution( idx_t npts, int part[], int part0 ) : impl_( new Implementation( npts, part, part0 ) ) {}

Distribution::~Distribution() {}

int Distribution::partition( const gidx_t gidx ) const {
    return impl_->partition( gidx );
}

const std::vector<int>& Distribution::partition() const {
    return impl_->partition();
}

idx_t Distribution::nb_partitions() const {
    return impl_->nb_partitions();
}

const int* Distribution::data() const {
    return impl_->data();
}

const std::vector<idx_t>& Distribution::nb_pts() const {
    return impl_->nb_pts();
}

idx_t Distribution::max_pts() const {
    return impl_->max_pts();
}

idx_t Distribution::min_pts() const {
    return impl_->min_pts();
}

const std::string& Distribution::type() const {
    return impl_->type();
}

std::ostream& operator<<( std::ostream& os, const Distribution& distribution ) {
    distribution.impl_->print( os );
    return os;
}

const Distribution::Implementation* Distribution::get() const {
    return impl_.get();
}

atlas::grid::Distribution::operator const std::vector<int>&() const {
    return *impl_;
}

}  // namespace grid
}  // namespace atlas
