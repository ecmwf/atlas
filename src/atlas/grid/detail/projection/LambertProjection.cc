#include <cmath>
#include "eckit/utils/MD5.h"
#include "atlas/grid/detail/projection/LambertProjection.h"
#include "atlas/util/Constants.h"

/*
Projection formula's for Lambert projection from "Map Projections: A Working Manual"

The origin of the xy-system is at (lon0,0)

*/

namespace {
  static double D2R(const double x) {
    return atlas::util::Constants::degreesToRadians()*x;
  }
  static double R2D(const double x) {
    return atlas::util::Constants::radiansToDegrees()*x;
  }
}

namespace atlas {
namespace grid {
namespace projection {

// constructors
LambertProjection::LambertProjection(const eckit::Parametrisation& params) {
  // check presence of radius
  if( ! params.get("radius",radius_) )
    radius_=util::Earth::radiusInMeters();
  // check presence of lat1 and lat2
  if( ! params.get("latitude1",lat1_) )
    throw eckit::BadParameter("latitude1 missing in Params",Here());
  if( ! params.get("latitude2",lat2_) )
    lat2_=lat1_;
  // check presence of lon0
  if( ! params.get("longitude0",lon0_) )
    throw eckit::BadParameter("longitude0 missing in Params",Here());

  setup();
}


void LambertProjection::setup() {
  // setup (derived) constants
  is_tangent_ = (lat1_==lat2_);
  if ( is_tangent_ ) {
    n_ = std::sin(D2R(lat1_));
  } else {
    n_ = std::log(std::cos(D2R(lat1_))/std::cos(D2R(lat2_)))
       / std::log(std::tan(D2R(45+lat2_*0.5))/std::tan(D2R(45.+lat1_*0.5)));
  }
  F_ = std::cos(D2R(lat1_))*std::pow(std::tan(D2R(45.+lat1_*0.5)),n_)/n_;
  rho0_  = radius_*F_;
  inv_n_ = 1./n_;
  sign_  = ( n_<0. ? -1. : 1. );
}

void LambertProjection::lonlat2xy(double crd[]) const {

  const double rho  = radius_*F_/std::pow(std::tan(D2R(45+crd[1]*0.5)),n_);
  const double theta = D2R( n_*(crd[0]-lon0_) );
  //eckit::geometry::reduceTo2Pi(theta); // bracket between 0 and 360
  // theta*=n_;
  crd[0] = rho*std::sin(theta);
  crd[1] = rho0_-rho*std::cos(theta);
}

// inverse projection
void LambertProjection::xy2lonlat(double crd[]) const {
  // auxiliaries
  const double y0    = rho0_-crd[1];
  const double rho   = sign_ * std::sqrt(crd[0]*crd[0]+y0*y0);
  const double theta = R2D(std::atan2(sign_*crd[0],sign_*y0));

  // longitude
  crd[0] = theta*inv_n_+lon0_;

  // latitude
  if (rho==0.) {
    crd[1] = sign_ * 90.;
  } else {
    crd[1] = 2.*R2D(std::atan(std::pow(radius_*F_/rho,inv_n_)))-90.;
  }
}

// specification
eckit::Properties LambertProjection::spec() const {
  eckit::Properties proj_spec;
  proj_spec.set("type",static_type());
  proj_spec.set("latitude1",lat1_);
  proj_spec.set("latitude2",lat2_);
  proj_spec.set("longitude0",lon0_);
  if( radius_ != util::Earth::radiusInMeters() ) proj_spec.set("radius",radius_);
  return proj_spec;
}

void LambertProjection::hash( eckit::MD5& md5 ) const {
  md5.add(static_type());
  md5.add(lat1_);
  md5.add(lat2_);
  md5.add(lon0_);
  md5.add(radius_);
}

register_BuilderT1(Projection,LambertProjection,LambertProjection::static_type());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

