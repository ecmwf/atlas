/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "DistributionImpl.h"

#include <algorithm>

#include "atlas/grid/detail/distribution/DistributionArray.h"


namespace atlas {
namespace grid {

DistributionImpl* atlas__GridDistribution__new( idx_t size, int part[], int part0 ) {
    return new detail::distribution::DistributionArray( 0, size, part, part0 );
}

void atlas__GridDistribution__delete( DistributionImpl* This ) {
    delete This;
}

void atlas__GridDistribution__nb_pts( DistributionImpl* This, idx_t nb_pts[] ) {
    const auto& nb_pts_ = This->nb_pts();
    std::copy( nb_pts_.begin(), nb_pts_.end(), &nb_pts[0] );
}

idx_t atlas__atlas__GridDistribution__nb_partitions( DistributionImpl* This ) {
    return This->nb_partitions();
}


}  // namespace grid
}  // namespace atlas
