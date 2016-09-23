#include "atlas/grid/projection/SchmidtProjection.h"

#include "atlas/util/Constants.h"
#include <cmath>

#define D2R(X) (util::Constants::degreesToRadians()*(X))
#define R2D(X) (util::Constants::radiansToDegrees()*(X))

namespace atlas {
namespace grid {
namespace projection {

// constructor
SchmidtProjection::SchmidtProjection(const eckit::Parametrisation& params) {
	if( ! params.get("stretching_factor",c_) )
    throw eckit::BadParameter("stretching_factor missing in Params",Here());
}

eckit::geometry::LLPoint2 SchmidtProjection::coords2lonlat(eckit::geometry::Point2 xy) {
	double lat=R2D(asin(cos(2*atan(1/c_*tan(acos(sin(D2R(xy[eckit::geometry::YY])))/2)))));
	
	return eckit::geometry::LLPoint2(xy[eckit::geometry::XX],lat);
}

eckit::geometry::Point2 SchmidtProjection::lonlat2coords(eckit::geometry::LLPoint2 ll) {

	double lat2=R2D(asin(cos(2*atan(c_*tan(acos(sin(D2R(ll.lat())))/2)))));

	return eckit::geometry::Point2(ll.lon(),lat2);
}

register_BuilderT1(Projection,SchmidtProjection,SchmidtProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

