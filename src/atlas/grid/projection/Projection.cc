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

Rotated::Rotated(const eckit::Parametrisation& p) {
  // get pole
  std::vector<double> pole(2);
  if( ! p.get("pole",pole) )
    throw eckit::BadParameter("pole missing in Params",Here());
  pole_.assign(pole[0],pole[1]);

  double latrp = D2R(90.0-pole_.lat());
  cos_latrp_ = std::cos(latrp);
  sin_latrp_ = std::sin(latrp);

}

Rotated::Rotated( const Rotated& rhs ) {
    pole_ = rhs.pole_;
    cos_latrp_ = rhs.cos_latrp_;
    sin_latrp_ = rhs.sin_latrp_;
}


void Rotated::rotate(eckit::geometry::LLPoint2 &P) const {
  // coordinates of the point P on a rotated sphere with specified pole

  double lon, lat, lonr, latr, lont, latt;
  double xt, yt, zt, x, y, z;
  double cos_lon, sin_lon, cos_lat, sin_lat;

  lon = D2R(P.lon());
  lat = D2R(P.lat());
  cos_lon = std::cos(lon);
  cos_lat = std::cos(lat);
  sin_lon = std::sin(lon);
  sin_lat = std::sin(lat);

  // cartesian coordinates
  x = cos_lon * cos_lat;
  y = sin_lon * cos_lat;
  z = sin_lat;

  // tilt
  xt = cos_latrp_*x + sin_latrp_*z;
  yt = y;
  zt = -sin_latrp_*x + cos_latrp_*z;

  // back to spherical coordinates
  lont=R2D(std::atan2(yt,xt));
  latt=R2D(std::asin(zt));

  // rotate
  lonr=lont+pole_.lon();
  latr=latt;

  P.assign(lonr,latr);
}

void Rotated::unrotate(eckit::geometry::LLPoint2 &P) const {
  // inverse operation of Projection::rotate

  double lon, lat, lont, latt;
  double xt, yt, zt, x, y, z;
  double cos_lont, sin_lont, cos_latt, sin_latt;

  // unrotate
  lont=D2R(P.lon()-pole_.lon());
  latt=D2R(P.lat());

  cos_lont  = std::cos(lont);
  cos_latt  = std::cos(latt);
  sin_lont  = std::sin(lont);
  sin_latt  = std::sin(latt);

  // cartesian coordinates
  xt = cos_lont * cos_latt;
  yt = sin_lont * cos_latt;
  zt = sin_latt;

  // untilt
  x = cos_latrp_*xt - sin_latrp_*zt;
  y = yt;
  z = sin_latrp_*xt + cos_latrp_*zt;

  // back to spherical coordinates
  lon=R2D(std::atan2(y,x));
  lat=R2D(std::asin(z));

  P.assign(lon,lat);
}

// specification
eckit::Properties Rotated::spec() const {
  eckit::Properties proj_spec;
  std::vector<double> p(2);
  p[0]=pole_.lon();
  p[1]=pole_.lat();
  proj_spec.set("projectionPole",eckit::makeVectorValue(p));
  return proj_spec;
}

}  // namespace projection
}  // namespace grid
}  // namespace atlas

