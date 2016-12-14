#include "atlas/grid/projection/RotatedLonLatProjection.h"

#include "atlas/util/Constants.h"
#include <cmath>

#define D2R(X) (util::Constants::degreesToRadians()*(X))
#define R2D(X) (util::Constants::radiansToDegrees()*(X))

namespace atlas {
namespace grid {
namespace projection {

// constructor
RotatedLonLatProjection::RotatedLonLatProjection(const eckit::Parametrisation& params) : LonLatProjection(params) {

	// get pole
  std::vector<double> p(2);
	if( ! params.get("pole",p) )
    throw eckit::BadParameter("pole missing in Params",Here());
  
  pole_.assign(p[0],p[1]);
  
}

eckit::geometry::LLPoint2 RotatedLonLatProjection::coords2lonlat(eckit::geometry::Point2 xy) {

	// point
	eckit::geometry::LLPoint2 P(xy[eckit::geometry::XX],xy[eckit::geometry::YY]);
	
	// perform rotation
	rotate_(P,pole_);

	return P;
}

eckit::geometry::Point2 RotatedLonLatProjection::lonlat2coords(eckit::geometry::LLPoint2 P) {

	// inverse rotation
	unrotate_(P,pole_);

	return eckit::geometry::Point2(P.lon(),P.lat());
}


// specification
eckit::Properties RotatedLonLatProjection::spec() const {
	eckit::Properties proj_spec;
	proj_spec.set("projectionType",virtual_projection_type_str());
	std::vector<double> p(2);
	p[0]=pole_.lon();
	p[1]=pole_.lat();
	proj_spec.set("projectionPole",eckit::makeVectorValue(p));
	return proj_spec;
}

register_BuilderT1(Projection,RotatedLonLatProjection,RotatedLonLatProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

