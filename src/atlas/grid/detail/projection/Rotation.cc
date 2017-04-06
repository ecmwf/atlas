#include <cmath>
#include "eckit/utils/MD5.h"
#include "atlas/grid/detail/projection/Rotation.h"
#include "atlas/util/Constants.h"

namespace atlas {
namespace grid {
namespace projection {

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
  std::vector<double> north_pole, south_pole;
  if( not p.get("north_pole",north_pole) ) {
    if( not p.get("south_pole",south_pole) ) {
      throw eckit::BadParameter("north_pole or south_pole missing in Params",Here());
    }
    north_pole.resize(2);
    north_pole[0] =  south_pole[0]-180.; // longitude
    north_pole[1] = -south_pole[1];      // latitude
  }
  pole_ = PointLonLat(north_pole[0],north_pole[1]);

  if( pole_.lon() == 0. && pole_.lat() == 90. ) rotated_ = false;

  double latrp = D2R(90.0-pole_.lat());
  cos_latrp_ = std::cos(latrp);
  sin_latrp_ = std::sin(latrp);
}

Rotated::Rotated( const Rotated& rhs ) {
    pole_ = rhs.pole_;
    cos_latrp_ = rhs.cos_latrp_;
    sin_latrp_ = rhs.sin_latrp_;
    rotated_ = rhs.rotated_;
}


void Rotated::rotate(double crd[]) const {
  // coordinates of the point P on a rotated sphere with specified pole
  double lon, lat, lont, latt;
  double xt, yt, zt, x, y, z;
  double cos_lon, sin_lon, cos_lat, sin_lat;

  lon = D2R(crd[0]);
  lat = D2R(crd[1]);
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
  crd[0]=lont+pole_.lon();
  crd[1]=latt;
}


void Rotated::unrotate(double crd[]) const {
  // inverse operation of Projection::rotate

  double lont, latt;
  double xt, yt, zt, x, y, z;
  double cos_lont, sin_lont, cos_latt, sin_latt;

  // unrotate
  lont=D2R(crd[0]-pole_.lon());
  latt=D2R(crd[1]);

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
  crd[0]=R2D(std::atan2(y,x));
  crd[1]=R2D(std::asin(z));
}

// specification
void Rotated::spec(eckit::Properties& s) const {
  std::vector<double> p(2);
  p[0]=pole_.lon();
  p[1]=pole_.lat();
  s.set("projectionPole",eckit::makeVectorValue(p));
}

void Rotated::hash( eckit::MD5& md5 ) const {
  md5.add("rotated");
  md5.add(pole_.lon());
  md5.add(pole_.lat());
}

void NotRotated::hash( eckit::MD5& md5 ) const {
  md5.add("not_rotated");
}


}  // namespace projection
}  // namespace grid
}  // namespace atlas

