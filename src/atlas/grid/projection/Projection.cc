#include "atlas/grid/projection/Projection.h"

#include "atlas/util/Constants.h"
#include <cmath>

#define D2R(X) (atlas::util::Constants::degreesToRadians()*(X))
#define R2D(X) (atlas::util::Constants::radiansToDegrees()*(X))

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
  return NULL;
}

void Projection::rotate_(eckit::geometry::LLPoint2 &P,const eckit::geometry::LLPoint2 &pole) const {
  // coordinates of the point P on a rotated sphere with specified pole

  double lon, lat, lonr, latr, lont, latt;
  double xt, yt, zt, x, y, z;

  lon=P.lon();
  lat=P.lat();

  // cartesian coordinates
  x=cos(D2R(lon))*cos(D2R(lat));
  y=sin(D2R(lon))*cos(D2R(lat));
  z=sin(D2R(lat));

  // tilt
  xt=cos(D2R(90.0-pole.lat()))*x + sin(D2R(90.0-pole.lat()))*z;
  yt=y;
  zt=-sin(D2R(90.0-pole.lat()))*x + cos(D2R(90.0-pole.lat()))*z;

  // back to spherical coordinates
  lont=R2D(atan2(yt,xt));
  latt=R2D(asin(zt));

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
  xt=cos(D2R(lont))*cos(D2R(latt));
  yt=sin(D2R(lont))*cos(D2R(latt));
  zt=sin(D2R(latt));

  // untilt
  x=cos(D2R(90.0-pole.lat()))*xt - sin(D2R(90.0-pole.lat()))*zt;
  y=yt;
  z=sin(D2R(90.0-pole.lat()))*xt + cos(D2R(90.0-pole.lat()))*zt;

  // back to spherical coordinates
  lon=R2D(atan2(y,x));
  lat=R2D(asin(z));

  P.assign(lon,lat);
}

}  // namespace projection
}  // namespace grid
}  // namespace atlas

