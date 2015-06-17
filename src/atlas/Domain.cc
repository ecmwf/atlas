/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/Domain.h"

#include "eckit/utils/MD5.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/BoundBox.h"

using eckit::FloatCompare;

namespace atlas {

Domain::Domain(const atlas::BoundBox& bbox) :
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
  double lat = p.lat();
  double lon = normalise(lon);

  return FloatCompare<double>::isGreaterApproxEqual(north_, lat) &&
         FloatCompare<double>::isGreaterApproxEqual(lat, south_) &&
         FloatCompare<double>::isGreaterApproxEqual(lon , west_) &&
         FloatCompare<double>::isGreaterApproxEqual(east_, lon);
}

Domain Domain::makeGlobal() { return Domain(90.,0.,-90.,360.); }

void Domain::normalise()
{
  while (east_ >= 360) {
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

}  // namespace atlas

