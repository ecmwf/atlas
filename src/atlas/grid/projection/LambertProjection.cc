#include <cmath>

#include "atlas/grid/projection/LambertProjection.h"
#include "atlas/util/Constants.h"

/*
Projection formula's for Lambert projection from "Map Projections: A Working Manual"

The origin of the xy-system is at (lon0,0)

*/

#define D2R(X) (util::Constants::degreesToRadians()*(X))
#define R2D(X) (util::Constants::radiansToDegrees()*(X))

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

// setup routine
void LambertProjection::setup() {
	// setup (derived) constants
	isTangent_=(lat1_==lat2_);
	if ( isTangent_ ) {
		n_=sin(D2R(lat1_));
	} else {
		n_=log(cos(D2R(lat1_))/cos(D2R(lat2_)))/log(tan(D2R(45+lat2_/2))/tan(D2R(45+lat1_/2)));
	}
	F_=cos(D2R(lat1_))*pow(tan(D2R(45+lat1_/2)),n_)/n_;
	rho0_=radius_*F_;
}

// projection
eckit::geometry::Point2 LambertProjection::lonlat2coords(eckit::geometry::LLPoint2 ll) const {

	double rho=radius_*F_/pow(tan(D2R(45+ll.lat()/2)),n_);
	double theta=ll.lon()-lon0_;
	eckit::geometry::reduceTo2Pi(theta);
	theta*=n_;
	double x=rho*sin(D2R(theta));
	double y=rho0_-rho*cos(D2R(theta));
	
	return eckit::geometry::Point2(x,y);
}

// inverse projection
eckit::geometry::LLPoint2 LambertProjection::coords2lonlat(eckit::geometry::Point2 xy) const {
	double x=xy[eckit::geometry::XX], y=xy[eckit::geometry::YY];
	
	// auxiliaries
	double rho=sqrt(x*x+(rho0_-y)*(rho0_-y));
	if (n_<0) rho=-rho;
	double theta;
	if (n_>0) {
		theta=R2D(atan2(x,rho0_-y));
	} else {
		theta=R2D(atan2(-x,y-rho0_));
	}
	
	// longitude
	double lon=theta/n_+lon0_;
	
	// latitude
	double lat;
	if (rho==0) {
		lat=( n_>0 ? 90 : -90 );
	} else {
		lat=2*R2D(atan(pow(radius_*F_/rho,1/n_)))-90;
	}

	return eckit::geometry::LLPoint2(lon,lat);
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

