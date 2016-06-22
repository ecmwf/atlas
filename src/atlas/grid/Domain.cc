/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grid/Domain.h"

#include "eckit/types/FloatCompare.h"
#include "eckit/utils/MD5.h"


namespace atlas {
namespace grid {


typedef eckit::FloatCompare<double> cmp;



void Domain::hash(eckit::MD5& md5) const {
    md5.add(north_);
    md5.add(west_);
    md5.add(south_);
    md5.add(east_);
}


bool Domain::isEmpty() const {
    return !cmp::isStrictlyGreater(north_,south_)
        || !cmp::isStrictlyGreater(east_,west_);
}


bool Domain::contains(double lon, double lat) const {
    // approximate comparisons include boundary coordinates
    lon = normalise(lon);
    return  cmp::isApproximatelyGreaterOrEqual(north_, lat) &&
            cmp::isApproximatelyGreaterOrEqual(lat, south_) &&
            cmp::isApproximatelyGreaterOrEqual(lon, west_)  &&
            cmp::isApproximatelyGreaterOrEqual(east_, lon);
}


void Domain::normalise() {
    ASSERT(north_ <= 90 && south_ >= -90);
    ASSERT(north_ >= south_);

    while (east_ >  360)  { east_ -= 360; west_ -= 360; }
    while (east_ < -180)  { east_ += 360; west_ += 360; }
    while (east_ < west_) { east_ += 360; }

    ASSERT(west_ <= east_);
}


double Domain::normalise(double lon) const {
    while (lon > east_) { lon -= 360; }
    while (lon < west_) { lon += 360; }
    return lon;
}


void Domain::print(std::ostream& os) const {
    os  << "Domain["
        <<  "N:" << north_
        << ",W:" << west_
        << ",S:" << south_
        << ",E:" << east_
        << "]";
}


} // namespace grid
} // namespace atlas

