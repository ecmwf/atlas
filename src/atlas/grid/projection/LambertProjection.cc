#include <cmath>
#include "atlas/grid/projection/LambertProjection.h"
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

// copy constructor
LambertProjection::LambertProjection( const LambertProjection& rhs ) : Projection( rhs ) {
  // copy fundamental data
  lat1_=rhs.lat1_;
  lat2_=rhs.lat2_;
  radius_=rhs.radius_;
  lon0_=rhs.lon0_;

  setup();
}

void LambertProjection::setup() {
  // setup (derived) constants
  isTangent_=(lat1_==lat2_);
  if ( isTangent_ ) {
    n_=std::sin(D2R(lat1_));
  } else {
    n_=std::log(std::cos(D2R(lat1_))/std::cos(D2R(lat2_)))/std::log(std::tan(D2R(45+lat2_*0.5))/std::tan(D2R(45.+lat1_*0.5)));
  }
  F_=std::cos(D2R(lat1_))*std::pow(std::tan(D2R(45.+lat1_*0.5)),n_)/n_;
  rho0_=radius_*F_;
  inv_n_ = 1./n_;
}

// clone method
Projection* LambertProjection::clone() const { return new LambertProjection(*this); }

// // projection
// eckit::geometry::Point2 LambertProjection::lonlat2coords(eckit::geometry::LLPoint2 ll) const {
//
//   double rho=radius_*F_/std::pow(std::tan(D2R(45+ll.lat()*0.5)),n_);
//   double theta=ll.lon()-lon0_;
//   eckit::geometry::reduceTo2Pi(theta);
//   theta*=n_;
//   double x=rho*std::sin(D2R(theta));
//   double y=rho0_-rho*std::cos(D2R(theta));
//
//   return eckit::geometry::Point2(x,y);
// }
//
// // inverse projection
// eckit::geometry::LLPoint2 LambertProjection::coords2lonlat(eckit::geometry::Point2 xy) const {
//   double x=xy[eckit::geometry::XX], y=xy[eckit::geometry::YY];
//
//   // auxiliaries
//   double rho=std::sqrt(x*x+(rho0_-y)*(rho0_-y));
//   if (n_<0.) rho = -rho;
//   double theta;
//   if (n_>0.) {
//     theta=R2D(std::atan2(x,rho0_-y));
//   } else {
//     theta=R2D(std::atan2(-x,y-rho0_));
//   }
//
//   // longitude
//   double lon=theta*inv_n_+lon0_;
//
//   // latitude
//   double lat;
//   if (rho==0.) {
//     lat=( n_>0. ? 90. : -90. );
//   } else {
//     lat=2.*R2D(std::atan(std::pow(radius_*F_/rho,inv_n_)))-90.;
//   }
//
//   return eckit::geometry::LLPoint2(lon,lat);
// }

void LambertProjection::lonlat2coords(double crd[]) const {

  double rho=radius_*F_/std::pow(std::tan(D2R(45+crd[1]*0.5)),n_);
  double theta=crd[0]-lon0_;
  eckit::geometry::reduceTo2Pi(theta); // bracket between 0 and 360
  theta*=n_;
  crd[0]=rho*std::sin(D2R(theta));
  crd[1]=rho0_-rho*std::cos(D2R(theta));
}

// inverse projection
void LambertProjection::coords2lonlat(double crd[]) const {
  // auxiliaries
  double rho=std::sqrt(crd[0]*crd[0]+(rho0_-crd[1])*(rho0_-crd[1]));
  if (n_<0.) rho = -rho;
  double theta;
  if (n_>0.) {
    theta=R2D(std::atan2(crd[0],rho0_-crd[1]));
  } else {
    theta=R2D(std::atan2(-crd[0],crd[1]-rho0_));
  }

  // longitude
  crd[0]=theta*inv_n_+lon0_;

  // latitude
  double lat;
  if (rho==0.) {
    crd[1]=( n_>0. ? 90. : -90. );
  } else {
    crd[1]=2.*R2D(std::atan(std::pow(radius_*F_/rho,inv_n_)))-90.;
  }
}

// specification
eckit::Properties LambertProjection::spec() const {
  eckit::Properties proj_spec;
  proj_spec.set("projectionType",virtual_projection_type_str());
  proj_spec.set("projectionLatitude1",lat1_);
  proj_spec.set("projectionLatitude2",lat2_);
  proj_spec.set("projectionLongitude0",lon0_);
  proj_spec.set("projectionRadius",radius_);
  return proj_spec;
}


register_BuilderT1(Projection,LambertProjection,LambertProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

