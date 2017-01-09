#include <cmath>

#include "atlas/grid/projection/RotatedMercatorProjection.h"
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
RotatedMercatorProjection::RotatedMercatorProjection(const eckit::Parametrisation& params) : MercatorProjection(params) {

	setup(params);

}

// copy constructor
RotatedMercatorProjection::RotatedMercatorProjection( const RotatedMercatorProjection& rhs ) : MercatorProjection(rhs) {
	pole_.assign(rhs.pole_[0],rhs.pole_[1]);
	lon0_=rhs.lon0_;
}

// clone method
RotatedMercatorProjection * RotatedMercatorProjection::clone() const  {
	return new RotatedMercatorProjection(*this);
}


void RotatedMercatorProjection::setup(const eckit::Parametrisation& params) {
  // check presence of pole
  std::vector<double> p(2);
	if( ! params.get("pole",p) )
    throw eckit::BadParameter("pole missing in Params",Here());
  pole_.assign(p[0],p[1]);
  
  MercatorProjection::setup(params);
  // rotation is in pole_
  lon0_=0.0;

}

// projection
eckit::geometry::Point2 RotatedMercatorProjection::lonlat2coords(eckit::geometry::LLPoint2 ll) const {

	// first unrotate
	eckit::geometry::LLPoint2 P(ll.lon(),ll.lat());
	unrotate_(P,pole_);
	
	// then project
	return MercatorProjection::lonlat2coords(P);
}

// inverse projection
eckit::geometry::LLPoint2 RotatedMercatorProjection::coords2lonlat(eckit::geometry::Point2 xy) const {

	// inverse projection
	eckit::geometry::LLPoint2 P=MercatorProjection::coords2lonlat(xy);
	
	// unrotate
	rotate_(P,pole_);
	
	// then project
	return P;
}

// specification
eckit::Properties RotatedMercatorProjection::spec() const {
	eckit::Properties proj_spec;
	proj_spec.set("projectionType",virtual_projection_type_str());
	proj_spec.set("projectionRadius",radius_);
	std::vector<double> p(2);
	p[0]=pole_.lon();
	p[1]=pole_.lat();
	proj_spec.set("projectionPole",eckit::makeVectorValue(p));
	return proj_spec;
}


register_BuilderT1(Projection,RotatedMercatorProjection,RotatedMercatorProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

