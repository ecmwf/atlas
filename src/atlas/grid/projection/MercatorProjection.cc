#include <cmath>

#include "atlas/grid/projection/MercatorProjection.h"
#include "atlas/util/Constants.h"

/*
Projection formula's for Mercator projection from "Map Projections: A Working Manual"

The origin of the xy-system is at (lon0,0)

*/

#define D2R(X) (util::Constants::degreesToRadians()*(X))
#define R2D(X) (util::Constants::radiansToDegrees()*(X))

namespace atlas {
namespace grid {
namespace projection {

// constructors
MercatorProjection::MercatorProjection(const eckit::Parametrisation& params) {
  setup(params);
}

// copy constructor
MercatorProjection::MercatorProjection( const MercatorProjection& rhs ) : Projection(rhs)  {
	radius_=rhs.radius_;
	lon0_=rhs.lon0_;
}

// clone method
MercatorProjection * MercatorProjection::clone() const  {
	return new MercatorProjection(*this);
}

// setup routine
void MercatorProjection::setup(const eckit::Parametrisation& params) {
	// check presence of radius
	if( ! params.get("radius",radius_) )
		radius_=util::Earth::radiusInMeters();
	// check presence of lon0
	if( ! params.get("longitude0",lon0_) ) 
    lon0_=0.0;
}

// projection
eckit::geometry::Point2 MercatorProjection::lonlat2coords(eckit::geometry::LLPoint2 ll) const {

	double x=radius_*(D2R(ll.lon()-lon0_));
	double y=radius_*log(tan(D2R(45+ll.lat()/2)));
	
	return eckit::geometry::Point2(x,y);
}

// inverse projection
eckit::geometry::LLPoint2 MercatorProjection::coords2lonlat(eckit::geometry::Point2 xy) const {
	double x=xy[eckit::geometry::XX], y=xy[eckit::geometry::YY];
	
	double lon=lon0_+R2D(x/radius_);
	double lat=2*R2D(atan(exp(y/radius_)))-90;

	return eckit::geometry::LLPoint2(lon,lat);
}

// specification
eckit::Properties MercatorProjection::spec() const {
	eckit::Properties proj_spec;
	proj_spec.set("projectionType",virtual_projection_type_str());
	proj_spec.set("projectionLongitude0",lon0_);
	proj_spec.set("projectionRadius",radius_);
	return proj_spec;
}

register_BuilderT1(Projection,MercatorProjection,MercatorProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

