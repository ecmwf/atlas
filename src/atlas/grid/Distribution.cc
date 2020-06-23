/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Distribution.h"

#include "eckit/utils/MD5.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/distribution/DistributionArray.h"
#include "atlas/grid/detail/distribution/SerialDistribution.h"

namespace atlas {
namespace grid {

using namespace detail::distribution;

Distribution::Distribution( const Grid& grid ) : Handle( new SerialDistribution{grid} ) {}

Distribution::Distribution( const Grid& grid, const Config& config ) :
    Handle( Partitioner( config ).partition( grid ).get() ) {}

Distribution::Distribution( const Grid& grid, const Partitioner& partitioner ) :
    Handle( partitioner.partition( grid ) ) {}

Distribution::Distribution( int nb_partitions, idx_t npts, int part[], int part0 ) :
    Handle( new DistributionArray( nb_partitions, npts, part, part0 ) ) {}

Distribution::Distribution( int nb_partitions, partition_t&& part ) :
    Handle( new DistributionArray( nb_partitions, std::move( part ) ) ) {}

Distribution::~Distribution() = default;

const std::vector<idx_t>& Distribution::nb_pts() const {
    return get()->nb_pts();
}

idx_t Distribution::max_pts() const {
    return get()->max_pts();
}

idx_t Distribution::min_pts() const {
    return get()->min_pts();
}

const std::string& Distribution::type() const {
    return get()->type();
}

std::ostream& operator<<( std::ostream& os, const Distribution& distribution ) {
    distribution.get()->print( os );
    return os;
}

void Distribution::hash( eckit::Hash& hash ) const {
    get()->hash( hash );
}

std::string Distribution::hash() const {
    eckit::MD5 h;
    hash( h );
    return h.digest();
}

}  // namespace grid
}  // namespace atlas
