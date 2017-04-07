#include <cmath>
#include "eckit/utils/MD5.h"
#include "atlas/grid/detail/projection/Rotation.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/runtime/Log.h"

#define LAM_IMPLEMENTATION 0

// Temporary option to  MIR is validated
#ifndef MIR_VALIDATE
#define MIR_VALIDATE 1
#endif

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
  
  
  const double theta = -(90.0 + spole_.lat()) * deg2rad;
  const double phi   = -spole_.lon()          * deg2rad;

  sin_theta = std::sin(theta);
  cos_theta = std::cos(theta);
  sin_phi   = std::sin(phi);
  cos_phi   = std::cos(phi);
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

    PointLonLat L(crd);
    L *= deg2rad;

    const double cos_lon = std::cos(L.lon());
    const double cos_lat = std::cos(L.lat());
    const double sin_lon = std::sin(L.lon());
    const double sin_lat = std::sin(L.lat());

    // cartesian coordinates
    const PointXYZ P(   cos_lon * cos_lat,
                        sin_lon * cos_lat,
                        sin_lat             );

    // Pt = Rot(z) * Rot(y) * P,   rotate about y axes then z
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

    PointXYZ Pt (  cos_theta*cos_phi*P.x()    + sin_phi*P.y()    + sin_theta*cos_phi*P.z() ,
                  -cos_theta*sin_phi*P.x()    + cos_phi*P.y()    - sin_theta*sin_phi*P.z() ,
                  -sin_theta        *P.x()                       + cos_theta        *P.z() );

    // Then convert back to 'normal' (lat,lon) using
    // Uses arc sin, to convert back to degrees, put in range -1 to 1 in case of slight rounding error
    // avoid error on calculating e.g. asin(1.00000001)
    if      (Pt.z() >  1.0)                 Pt.z()  =  1.0;
    else if (Pt.z() < -1.0)                 Pt.z()  = -1.0;
    if      (std::abs(Pt.y()) <  1.e-15 )   Pt.y() *=  0.0;
    
    PointLonLat Lt (  std::atan2(Pt.y(), Pt.x()) , 
                      std::asin (Pt.z())         );
    Lt *= rad2deg;

#ifdef MIR_VALIDATE
    // Still get a very small rounding error, round to 6 decimal places
    Lt.lat() = roundf( Lt.lat() * 1000000.0 )/1000000.0;
    Lt.lon() = roundf( Lt.lon() * 1000000.0 )/1000000.0;
#endif
    
    Lt.lon() -= angle_;

    // Make sure longitude is in range
    //     while ( Lt.lon() <  lonmin_) Lt.lon() += 360.0;
    //     while ( Lt.lon() >= lonmax_) Lt.lon() -= 360.0;

    crd[LON] = Lt.lon();
    crd[LAT] = Lt.lat();
}


void Rotated::unrotate(double crd[]) const {

    // See: http://rbrundritt.wordpress.com/2008/10/14/conversion-between-spherical-and-cartesian-coordinates-systems/
    // First convert the data point from spherical lat lon to (x',y',z') using:
    PointLonLat Lt (  crd[LON] + angle_ ,
                      crd[LAT]          );
    Lt *= deg2rad;
    
    // Cartesian coordinates
    const double cos_lont = std::cos(Lt.lon());
    const double cos_latt = std::cos(Lt.lat());
    const double sin_lont = std::sin(Lt.lon());
    const double sin_latt = std::sin(Lt.lat());
    
    const PointXYZ Pt ( cos_lont * cos_latt ,
                        sin_lont * cos_latt ,
                        sin_latt            );
    
    // Assume right hand rule, rotate about z axes and then y
    // P = Rot(y) * Rot(z) * Pt
    // x   (  cos(ϑ), 0, -sin(ϑ)) ( cos(φ), -sin(φ), 0) (xt)
    // y = (  0     , 1,  0     ) ( sin(φ), cos(φ),  0) (yt)
    // z   ( sin(ϑ), 0,   cos(ϑ)) ( 0     , 0     ,  1) (zt)

    // Expanded
    // x   ( cos(ϑ)cos(φ) , -cos(ϑ)sin(φ) , -sin(ϑ)) (xt)
    // y = ( sin(φ)       ,  cos(φ)       ,  0     ).(yt)
    // z   ( sin(ϑ) cos(φ), -sin(ϑ) sin(φ),  cos(ϑ)) (zt)

    PointXYZ P (
        cos_theta*cos_phi*Pt.x()   - cos_theta*sin_phi*Pt.y()    - sin_theta*Pt.z() ,
                  sin_phi*Pt.x()   +           cos_phi*Pt.y()                       ,
        sin_theta*cos_phi*Pt.x()   - sin_theta*sin_phi*Pt.y()    + cos_theta*Pt.z() );

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
    if      (P.z() >  1.0)                P.z() =  1.0;
    else if (P.z() < -1.0)                P.z() = -1.0;
    if      (std::abs(P.y()) <  1.e-15 )  P.y() *= 0.;

    PointLonLat L( std::atan2(P.y(), P.x()) ,
                    std::asin(P.z())        );
    L *= rad2deg;

#ifdef MIR_VALIDATE
    // Still get a very small rounding error, round to 6 decimal places
    L.lat() = roundf( L.lat() * 1000000.0 )/1000000.0;
    L.lon() = roundf( L.lon() * 1000000.0 )/1000000.0;
#endif

    // Make sure ret_lon is in range
    //     while (L.lon() <  lonmin_) L.lon() += 360.0;
    //     while (L.lon() >= lonmax_) L.lon() -= 360.0;

    crd[0] = L.lon();
    crd[1] = L.lat();
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

