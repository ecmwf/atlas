/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "DistributionFunction.h"

#include <algorithm>
#include <ostream>

#include "eckit/types/Types.h"
#include "eckit/utils/Hash.h"

#include "atlas/grid/Grid.h"

namespace atlas {
namespace grid {
namespace detail {
namespace distribution {

void DistributionFunction::print( std::ostream& s ) const {
    auto print_partition = [&]( std::ostream& s ) {
        eckit::output_list<int> list_printer( s );
        for ( gidx_t i = 0; i < size_; i++ ) {
            list_printer.push_back( partition( i ) );
        }
    };
    s << "Distribution( "
      << "type: " << type_ << ", nb_points: " << size_ << ", nb_partitions: " << nb_pts_.size() << ", parts : ";
    print_partition( s );
}

void DistributionFunction::hash( eckit::Hash& hash ) const {
    for ( gidx_t i = 0; i < size_; i++ ) {
        hash.add( partition( i ) );
    }
}

}  // namespace distribution
}  // namespace detail
}  // namespace grid
}  // namespace atlas
