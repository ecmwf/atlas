#include "atlas/grid/projections/LonLatProjection.h"

namespace atlas {
namespace grid {

LonLatProjection::LonLatProjection(const eckit::Parametrisation& params) {
}

eckit::geometry::LLPoint2 LonLatProjection::coords2lonlat(eckit::geometry::Point2 xy) {
	
	return eckit::geometry::LLPoint2(xy[eckit::geometry::XX],xy[eckit::geometry::YY]);
}

eckit::geometry::Point2 LonLatProjection::lonlat2coords(eckit::geometry::LLPoint2 ll) {

	return eckit::geometry::Point2(ll.lon(),ll.lat());
}

register_BuilderT1(Projection,LonLatProjection,LonLatProjection::className());

}  // namespace grid
}  // namespace atlas

