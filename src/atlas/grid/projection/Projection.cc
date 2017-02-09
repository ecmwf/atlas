#include <cmath>

#include "atlas/grid/projection/Projection.h"
#include "atlas/util/Constants.h"

namespace atlas {
namespace grid {
namespace projection {

Projection* Projection::create() {
  // default: no projection, i.e. stay in (lon,lat)-space
  util::Config projParams;
  projParams.set("projectionType","atlas.LonLatProjection");
  return Projection::create(projParams);
}

Projection* Projection::create(const eckit::Parametrisation& p) {
  std::string projectionType;
  if (p.get("projectionType",projectionType)) {
    return eckit::Factory<Projection>::instance().get(projectionType).create(p);
  }

  // should return error here
  throw eckit::BadParameter("projectionType missing in Params",Here());
}



namespace {
  static double D2R(const double x) {
    return atlas::util::Constants::degreesToRadians()*x;
  }
  static double R2D(const double x) {
    return atlas::util::Constants::radiansToDegrees()*x;
  }
}

void Projection::rotate_(eckit::geometry::LLPoint2 &P,const eckit::geometry::LLPoint2 &pole) const {
  // coordinates of the point P on a rotated sphere with specified pole

  double lon, lat, lonr, latr, lont, latt;
  double xt, yt, zt, x, y, z;

  lon=P.lon();
  lat=P.lat();

  // cartesian coordinates
  x=std::cos(D2R(lon))*std::cos(D2R(lat));
  y=std::sin(D2R(lon))*std::cos(D2R(lat));
  z=std::sin(D2R(lat));

  // tilt
  xt=std::cos(D2R(90.0-pole.lat()))*x + std::sin(D2R(90.0-pole.lat()))*z;
  yt=y;
  zt=-std::sin(D2R(90.0-pole.lat()))*x + std::cos(D2R(90.0-pole.lat()))*z;

  // back to spherical coordinates
  lont=R2D(std::atan2(yt,xt));
  latt=R2D(std::asin(zt));

  // rotate
  lonr=lont+pole.lon();
  latr=latt;

  P.assign(lonr,latr);
}

void Projection::unrotate_(eckit::geometry::LLPoint2 &P,const eckit::geometry::LLPoint2 &pole) const {
  // inverse operation of Projection::rotate

  double lon, lat, lonr, latr, lont, latt;
  double xt, yt, zt, x, y, z;

  lonr=P.lon();
  latr=P.lat();

  // unrotate
  lont=lonr-pole.lon();
  latt=latr;

  // cartesian coordinates
  xt=std::cos(D2R(lont))*std::cos(D2R(latt));
  yt=std::sin(D2R(lont))*std::cos(D2R(latt));
  zt=std::sin(D2R(latt));

  // untilt
  x=std::cos(D2R(90.0-pole.lat()))*xt - std::sin(D2R(90.0-pole.lat()))*zt;
  y=yt;
  z=std::sin(D2R(90.0-pole.lat()))*xt + std::cos(D2R(90.0-pole.lat()))*zt;

  // back to spherical coordinates
  lon=R2D(std::atan2(y,x));
  lat=R2D(std::asin(z));

  P.assign(lon,lat);
}

}  // namespace projection
}  // namespace grid
}  // namespace atlas

