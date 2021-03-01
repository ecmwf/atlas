/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "SerialDistribution.h"

#include <algorithm>
#include <ostream>

#include "atlas/grid/Grid.h"

namespace atlas {
namespace grid {
namespace detail {
namespace distribution {

SerialDistribution::SerialDistribution( const Grid& grid ) : DistributionFunctionT<SerialDistribution>( grid ) {
    type_          = "serial";
    nb_partitions_ = 1;
    size_          = grid.size();
    nb_pts_.resize( nb_partitions_, grid.size() );
    max_pts_ = *std::max_element( nb_pts_.begin(), nb_pts_.end() );
    min_pts_ = *std::min_element( nb_pts_.begin(), nb_pts_.end() );
}

}  // namespace distribution
}  // namespace detail
}  // namespace grid
}  // namespace atlas
