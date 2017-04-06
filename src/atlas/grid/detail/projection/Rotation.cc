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
#if 0
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
#else
    const double
            degree_to_radian_ = D2R(1.),
            radian_to_degree_ = R2D(1.),
            south_pole_lon =  pole_[0] + 180.,
            south_pole_lat = -pole_[1],
            south_pole_rot_angle = 0,
            lonmin_ = -180,
            lonmax_ =  180;

    // See: http://gis.stackexchange.com/questions/10808/lon-lat-transformation/14445
    // First convert the data point from spherical lat lon to (x',y',z') using:
    double latr = crd[1] * degree_to_radian_ ;
    double lonr = crd[0] * degree_to_radian_ ;
    double xd = cos(lonr)*cos(latr);
    double yd = sin(lonr)*cos(latr);
    double zd = sin(latr);

    // P' = Rot(z) * Rot(y) * Pv,   rotate about y axes then Z
    // Since we're undoing the rotation described in the definition of the coordinate system,
    // we first rotate by ϑ = -(90 + south_pole_lat) around the y' axis (along the rotated Greenwich meridian)
    // and then by φ = -south_pole_lon = +15 degrees around the z axis):
    // x   ( cos(φ), sin(φ), 0) (  cos(ϑ), 0, sin(ϑ)) (x')
    // y = (-sin(φ), cos(φ), 0).(  0     , 1, 0     ).(y')
    // z   ( 0     , 0     , 1) ( -sin(ϑ), 0, cos(ϑ)) (z')

    // Expanded
    // x =  cos(ϑ) cos(φ) x' + sin(φ) y' + sin(ϑ) cos(φ) z'
    // y = -cos(ϑ) sin(φ) x' + cos(φ) y' - sin(ϑ) sin(φ) z'
    // z = -sin(ϑ) x' + cos(ϑ) z'

    double t = -(90.0 + south_pole_lat);
    double o = -south_pole_lon;

    double sin_t = sin(degree_to_radian_ * t);
    double cos_t = cos(degree_to_radian_ * t);
    double sin_o = sin(degree_to_radian_ * o);
    double cos_o = cos(degree_to_radian_ * o);

    double x = cos_t*cos_o*xd + sin_o*yd + sin_t*cos_o*zd;
    double y = -cos_t*sin_o*xd + cos_o*yd - sin_t*sin_o*zd;
    double z = -sin_t*xd + cos_t*zd;

    // Then convert back to 'normal' (lat,lon) using
    // Uses arc sin, to convert back to degrees, put in range -1 to 1 in case of slight rounding error
    // avoid error on calculating e.g. asin(1.00000001)
    if (z > 1.0)  z = 1.0;
    if (z < -1.0) z = -1.0;

    double ret_lat = asin(z)*radian_to_degree_;
    double ret_lon = atan2(y, x)*radian_to_degree_;

//    // Still get a very small rounding error, round to 6 decimal places
//    ret_lat = roundf( ret_lat * 1000000.0 )/1000000.0;
//    ret_lon = roundf( ret_lon * 1000000.0 )/1000000.0;

    ret_lon -= south_pole_rot_angle;

    // Make sure ret_lon is in range
    while (ret_lon < lonmin_) ret_lon += 360.0;
    while (ret_lon >= lonmax_) ret_lon -= 360.0;

    crd[0] = ret_lon;
    crd[1] = ret_lat;
#endif
}


void Rotated::unrotate(double crd[]) const {
#if 0
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
#else
    const double
            degree_to_radian_ = D2R(1.),
            radian_to_degree_ = R2D(1.),
            south_pole_lon =  pole_[0] + 180.,
            south_pole_lat = -pole_[1],
            south_pole_rot_angle = 0,
            lonmin_ = -180,
            lonmax_ =  180;

    // See: http://rbrundritt.wordpress.com/2008/10/14/conversion-between-spherical-and-cartesian-coordinates-systems/
    // First convert the data point from spherical lat lon to (x',y',z') using:
    double latr = crd[1] * degree_to_radian_ ;
    double lonr = crd[0] * degree_to_radian_ ;
    double xd = cos(lonr)*cos(latr);
    double yd = sin(lonr)*cos(latr);
    double zd = sin(latr);

    // Assume right hand rule, rotate about z axes and then y
    // P' = Rot(y) * Rot(z) * Pv
    // x   (  cos(ϑ), 0, -sin(ϑ)) ( cos(φ), -sin(φ), 0) (x')
    // y = (  0     , 1,  0     ) ( sin(φ), cos(φ),  0) (y')
    // z   ( sin(ϑ), 0,   cos(ϑ)) ( 0     , 0     ,  1) (z')

    // Expanded
    // x   ( cos(ϑ)cos(φ) , -cos(ϑ)sin(φ) , -sin(ϑ)) (x')
    // y = ( sin(φ)       ,  cos(φ)       ,  0     ).(y')
    // z   ( sin(ϑ) cos(φ), -sin(ϑ) sin(φ),  cos(ϑ)) (z')

    double t = -(90.0 + south_pole_lat);
    double o = -south_pole_lon + south_pole_rot_angle;

    double sin_t = sin(degree_to_radian_ * t);
    double cos_t = cos(degree_to_radian_ * t);
    double sin_o = sin(degree_to_radian_ * o);
    double cos_o = cos(degree_to_radian_ * o);

    double x = cos_t*cos_o*xd - cos_t*sin_o*yd - sin_t*zd;
    double y = sin_o*xd + cos_o*yd  ;
    double z = sin_t*cos_o*xd - sin_t*sin_o*yd + cos_t*zd;

    // Then convert back to 'normal' (lat,lon)
    // z = r.cosϑ  ( r is earths radius, assume 1)
    // r = sqrt(x.x + y.y + z.z)
    // z = r.cosϑ  => ϑ = cos-1(z/r)
    // ϑ = con-1( z/ sqrt(x.x + y.y + z.z) )
    // By rearranging the formulas for the x and y components we can solve the value for the Φ angle.
    //    y                  x
    //   ---     = sinφ =   ----
    //   r.sinϑ             r.cosϑ
    //
    //  y/x = sinϑ/cosϑ = tanϑ
    //
    //  ϑ = tan-1(y/x)   => lon = atan2(y, x)


    // Uses arc sin, to convert back to degrees, put in range -1 to 1 in case of slight rounding error
    // avoid error on calculating e.g. asin(1.00000001)
    if (z > 1.0)  z = 1.0;
    if (z < -1.0) z = -1.0;
    double ret_lat = asin(z) * radian_to_degree_;
    double ret_lon = atan2(y, x) * radian_to_degree_;

//    // Still get a very small rounding error, round to 6 decimal places
//    ret_lat = roundf( ret_lat * 1000000.0 )/1000000.0;
//    ret_lon = roundf( ret_lon * 1000000.0 )/1000000.0;

    // Make sure ret_lon is in range
    while (ret_lon < lonmin_) ret_lon += 360.0;
    while (ret_lon >= lonmax_) ret_lon -= 360.0;

    crd[0] = ret_lon;
    crd[1] = ret_lat;
#endif
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

