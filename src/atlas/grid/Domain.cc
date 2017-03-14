/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/Domain.h"

#include <algorithm>
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/MD5.h"


namespace atlas {
namespace grid {

//----------------------------------------------------------------------------------------------------------------------

using eckit::types::is_approximately_equal;
using eckit::types::is_strictly_greater;
using eckit::types::is_approximately_greater_or_equal;


void Domain::hash(eckit::MD5& md5) const {
    md5.add(north_);
    md5.add(west_);
    md5.add(south_);
    md5.add(east_);
}


bool Domain::contains(double lon, double lat) const {
    // approximate comparisons include boundary coordinates
    lon = normalise(lon);
    return  is_approximately_greater_or_equal(north_, lat) &&
            is_approximately_greater_or_equal(lat, south_) &&
            is_approximately_greater_or_equal(lon, west_)  &&
            is_approximately_greater_or_equal(east_, lon);
}


#if 0
Domain Domain::intersect(const Domain& other) const {
    return !intersects(other)? makeEmpty()
           : Domain(
               std::min(n, other.n),
               std::min(w, other.w),
               std::max(s, other.s),
               std::max(e, other.e) );
}


bool Domain::intersects(const Domain& other) const {
    ASSERT(other.n >= other.s);
    ASSERT(other.w <= other.e);

    // strict comparisons ensure resulting areas have 2-dimensionality
    return  is_strictly_greater(e, other.w) &&
            is_strictly_greater(other.e, w) &&
            is_strictly_greater(s, other.n) &&
            is_strictly_greater(other.s, n);
}
#endif


void Domain::normalise() {
    ASSERT(north_ <= 90 && south_ >= -90);
    ASSERT(north_ >= south_);

    while (east_ >  360) {
        east_ -= 360;
        west_ -= 360;
    }
    while (east_ < -180) {
        east_ += 360;
        west_ += 360;
    }
    while (east_ < west_) {
        east_ += 360;
    }

    ASSERT(west_ <= east_);
}


double Domain::normalise(double lon) const {

    // check if lon is coincident with the West meridian
    double increment = 360;
    if (is_approximately_equal(long((lon - west_)/increment) * increment, lon - west_)) {
        return west_;
    }

    // check if lon is in ]West, East[
    while (is_strictly_greater(west_, lon)) {
        lon += 360;
    }
    while (is_strictly_greater(lon, east_)) {
        lon -= 360;
    }
    return lon;
}


void Domain::print(std::ostream& os) const {
    os  << "Domain["
        <<  "N:" << north_
        << ",W:" << west_
        << ",S:" << south_
        << ",E:" << east_
        << ",isGlobal=" << isGlobal()
        << "]";
}

//----------------------------------------------------------------------------------------------------------------------


} // namespace grid
} // namespace atlas
