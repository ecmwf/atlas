/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/utils/MD5.h"
#include "eckit/types/FloatCompare.h"
#include "atlas/grid/Domain.h"
#include "atlas/grid/BoundBox.h"

using eckit::FloatCompare;

namespace atlas {
namespace grid {

Domain::Domain(const atlas::grid::BoundBox& bbox) :
    north_(bbox.north()),
    west_(bbox.west()),
    south_(bbox.south()),
    east_(bbox.east())
{
    normalise();
}

Domain::Domain(double north, double west, double south, double east) :
  north_(north),
  west_(west),
  south_(south),
  east_(east)
{
  normalise();
}

Domain::~Domain()
{
}

void Domain::hash(eckit::MD5& md5) const
{
  md5.add(north_);
  md5.add(west_);
  md5.add(south_);
  md5.add(east_);
}

bool Domain::contains( const eckit::geometry::LLPoint2& p ) const
{
    return contains(p.lon(),p.lat());
}

bool Domain::contains(double lon, double lat) const
{
    lon = normalise(lon);

    return FloatCompare<double>::isApproximatelyGreaterOrEqual(north_, lat) &&
           FloatCompare<double>::isApproximatelyGreaterOrEqual(lat, south_) &&
           FloatCompare<double>::isApproximatelyGreaterOrEqual(lon , west_) &&
           FloatCompare<double>::isApproximatelyGreaterOrEqual(east_, lon);
}

Domain Domain::makeGlobal() {
    return Domain(90.,0.,-90.,360.);
}

bool Domain::global() const
{
    // This logic should not be changed
    // We should not use increments to test for globalness

    return north_ == 90. && south_ == -90. && (east_ - west_) == 360.;
}

void Domain::normalise()
{
  while (east_ > 360) {
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

  ASSERT(north_ <= 90 && south_ <= 90 && north_ >= -90 && south_ >= -90);
  ASSERT(north_ >= south_);
  ASSERT(west_ <= east_);
}

double Domain::normalise(double lon) const
{
  while (lon > east_) {
      lon -= 360;
  }

  while (lon < west_) {
      lon += 360;
  }
  return lon;
}

void Domain::print(std::ostream& os) const
{
  os << "Domain("
     << "N:" << north_
     << ",W:" << west_
     << ",S:" << south_
     << ",E:" << east_
     << ")";
}

} // namespace grid
} // namespace atlas

