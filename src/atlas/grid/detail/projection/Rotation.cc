#include <cmath>
#include "eckit/utils/MD5.h"
#include "atlas/grid/detail/projection/Rotation.h"
#include "atlas/util/Constants.h"
#include "atlas/runtime/Log.h"

#define LAM_IMPLEMENTATION 0


namespace atlas {
namespace grid {
namespace projection {

namespace {
  static double deg2rad = atlas::util::Constants::degreesToRadians();
  static double rad2deg = atlas::util::Constants::radiansToDegrees();
}

void Rotated::print( std::ostream& out ) const {
  out << "north_pole:"<<npole_<<", south_pole:"<<spole_<<", rotation_angle:"<<angle_; 
}
std::ostream& operator<< (std::ostream& out, const Rotated& r) {
  r.print(out);
  return out;
}

Rotated::Rotated(const eckit::Parametrisation& p) {
  // get pole
  std::vector<double> pole(2);
  if( p.get("north_pole",pole) ) {
    npole_ = PointLonLat(pole.data());
    spole_ = PointLonLat(npole_.lon()+180.,npole_.lat()-180.);
    if( spole_.lat() < -90 ) spole_.lon() -= 180.;
  } else if( p.get("south_pole",pole) ) {
    spole_ = PointLonLat(pole.data());
    npole_ = PointLonLat(spole_.lon()-180.,spole_.lat()+180.);
    if( npole_.lat() > 90 ) npole_.lon() += 180.;
  } else {
    throw eckit::BadParameter("north_pole or south_pole missing in Params",Here());
  }

  // TODO!!!!
  if( npole_.lon() == 0. && npole_.lat() == 90. ) rotated_ = false;

  double latrp = (90.0-npole_.lat()) * deg2rad;
  cos_latrp_ = std::cos(latrp);
  sin_latrp_ = std::sin(latrp);
  
  angle_ = 0.;
  p.get("rotation_angle",angle_);
}

Rotated::Rotated( const Rotated& rhs ) {
    npole_ = rhs.npole_;
    spole_ = rhs.spole_;
    cos_latrp_ = rhs.cos_latrp_;
    sin_latrp_ = rhs.sin_latrp_;
    rotated_ = rhs.rotated_;
    angle_ = rhs.angle_;
    lonmin_ = rhs.lonmin_;
}

#if OLD_IMPLEMENTATION

void Rotated::rotate(double crd[]) const {
  // coordinates of the point P on a rotated sphere with specified pole
  double lon, lat, lont, latt;
  double xt, yt, zt, x, y, z;
  double cos_lon, sin_lon, cos_lat, sin_lat;

  lon = crd[0] * deg2rad;
  lat = crd[1] * deg2rad;
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
  lont=std::atan2(yt,xt)) * rad2deg;
  latt=std::asin(zt))     * rad2deg;

  // rotate
  crd[0]=lont+npole_.lon();
  crd[1]=latt;
}


void Rotated::unrotate(double crd[]) const {
  // inverse operation of Projection::rotate

  double lont, latt;
  double xt, yt, zt, x, y, z;
  double cos_lont, sin_lont, cos_latt, sin_latt;

  // unrotate
  lont=(crd[0]-npole_.lon()) * deg2rad;
  latt=(crd[1])              * deg2rad;

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
  crd[0]= std::atan2(y,x) * rad2deg;
  crd[1]= std::asin(z)    * rad2deg;
}

#else

void Rotated::rotate(double crd[]) const {

    // See: http://gis.stackexchange.com/questions/10808/lon-lat-transformation/14445
    // First convert the data point from spherical lon lat to (x',y',z') using:

    double lon = crd[0] * deg2rad;
    double lat = crd[1] * deg2rad;

    const double cos_lon = std::cos(lon);
    const double cos_lat = std::cos(lat);
    const double sin_lon = std::sin(lon);
    const double sin_lat = std::sin(lat);

    // cartesian coordinates
    const double x = cos_lon * cos_lat;
    const double y = sin_lon * cos_lat;
    const double z = sin_lat;

    // P' = Rot(z) * Rot(y) * Pv,   rotate about y axes then z
    // Since we're undoing the rotation described in the definition of the coordinate system,
    // we first rotate by ϑ = -(90 + spole_.lat()) around the y axis (along the rotated Greenwich meridian)
    // and then by φ = -spole_.lon() degrees around the z axis):
    // xt   ( cos(φ), sin(φ), 0) (  cos(ϑ), 0, sin(ϑ)) (x)
    // yt = (-sin(φ), cos(φ), 0).(  0     , 1, 0     ).(y)
    // zt   ( 0     , 0     , 1) ( -sin(ϑ), 0, cos(ϑ)) (z)

    // Expanded
    // xt =  cos(ϑ) cos(φ) x + sin(φ) y + sin(ϑ) cos(φ) z
    // yt = -cos(ϑ) sin(φ) x + cos(φ) y - sin(ϑ) sin(φ) z
    // zt = -sin(ϑ)        x            + cos(ϑ)        z

    const double theta = -(90.0 + spole_.lat()) * deg2rad;
    const double phi   = -spole_.lon()          * deg2rad;

    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);
    const double sin_phi   = std::sin(phi);
    const double cos_phi   = std::cos(phi);

    double xt =  cos_theta*cos_phi*x    + sin_phi*y    + sin_theta*cos_phi*z;
    double yt = -cos_theta*sin_phi*x    + cos_phi*y    - sin_theta*sin_phi*z;
    double zt = -sin_theta        *x                   + cos_theta        *z;

    // Then convert back to 'normal' (lat,lon) using
    // Uses arc sin, to convert back to degrees, put in range -1 to 1 in case of slight rounding error
    // avoid error on calculating e.g. asin(1.00000001)
    if      (zt >  1.0)   zt =  1.0;
    else if (zt < -1.0)   zt = -1.0;
    if      (std::abs(yt) <  1.e-15 ) yt *= 0.;
    

    double lont = std::atan2(yt, xt) * rad2deg;
    double latt = std::asin(zt)      * rad2deg;

    // Still get a very small rounding error, round to 6 decimal places
    //    latt = roundf( latt * 1000000.0 )/1000000.0;
    //    lont = roundf( lont * 1000000.0 )/1000000.0;

    lont -= angle_;

    // Make sure ret_lon is in range
    while ( lont <  lonmin_) lont += 360.0;
    while ( lont >= lonmax_) lont -= 360.0;

    crd[0] = lont;
    crd[1] = latt;
}


void Rotated::unrotate(double crd[]) const {

    // See: http://rbrundritt.wordpress.com/2008/10/14/conversion-between-spherical-and-cartesian-coordinates-systems/
    // First convert the data point from spherical lat lon to (x',y',z') using:
    const double lont = (crd[0] + angle_) * deg2rad;
    const double latt = crd[1]            * deg2rad;
    
    // Cartesian coordinates
    const double cos_lont = std::cos(lont);
    const double cos_latt = std::cos(latt);
    const double sin_lont = std::sin(lont);
    const double sin_latt = std::sin(latt);
    const double xt = cos_lont * cos_latt;
    const double yt = sin_lont * cos_latt;
    const double zt = sin_latt;
    
    // Assume right hand rule, rotate about z axes and then y
    // P' = Rot(y) * Rot(z) * Pv
    // x   (  cos(ϑ), 0, -sin(ϑ)) ( cos(φ), -sin(φ), 0) (xt)
    // y = (  0     , 1,  0     ) ( sin(φ), cos(φ),  0) (yt)
    // z   ( sin(ϑ), 0,   cos(ϑ)) ( 0     , 0     ,  1) (zt)

    // Expanded
    // x   ( cos(ϑ)cos(φ) , -cos(ϑ)sin(φ) , -sin(ϑ)) (xt)
    // y = ( sin(φ)       ,  cos(φ)       ,  0     ).(yt)
    // z   ( sin(ϑ) cos(φ), -sin(ϑ) sin(φ),  cos(ϑ)) (zt)

    const double theta = -(90.0 + spole_.lat()) * deg2rad;
    const double phi   = -spole_.lon()          * deg2rad;

    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);
    const double sin_phi   = std::sin(phi);
    const double cos_phi   = std::cos(phi);

    double x = cos_theta*cos_phi*xt   - cos_theta*sin_phi*yt    - sin_theta*zt;
    double y =           sin_phi*xt   +           cos_phi*yt;
    double z = sin_theta*cos_phi*xt   - sin_theta*sin_phi*yt    + cos_theta*zt;

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
    if      (z >  1.0) z =  1.0;
    else if (z < -1.0) z = -1.0;
    if      (std::abs(y) <  1.e-15 ) y *= 0.;

    double lon = std::atan2(y, x) * rad2deg;
    double lat = std::asin(z)     * rad2deg;

    // Still get a very small rounding error, round to 6 decimal places
    //    ret_lat = roundf( ret_lat * 1000000.0 )/1000000.0;
    //    ret_lon = roundf( ret_lon * 1000000.0 )/1000000.0;

    // Make sure ret_lon is in range
    // while (lon <  lonmin_) lon += 360.0;
    // while (lon >= lonmax_) lon -= 360.0;
    
    crd[0] = lon;
    crd[1] = lat;
}

#endif

// specification
void Rotated::spec(eckit::Properties& s) const {
  std::vector<double> npole{ npole_.lon(), npole_.lat() };
  std::vector<double> spole{ spole_.lon(), spole_.lat() };
  s.set("north_pole",eckit::makeVectorValue(npole));
  s.set("south_pole",eckit::makeVectorValue(spole));
  s.set("rotation_angle",angle_);
}

void Rotated::hash( eckit::MD5& md5 ) const {
  md5.add("rotated");
  md5.add(spole_.lon());
  md5.add(spole_.lat());
  md5.add(angle_);
  md5.add(lonmin_);
}

void NotRotated::hash( eckit::MD5& md5 ) const {
  md5.add("not_rotated");
}


}  // namespace projection
}  // namespace grid
}  // namespace atlas

