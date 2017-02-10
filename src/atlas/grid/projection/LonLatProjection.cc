#include "atlas/grid/projection/LonLatProjection.h"

namespace atlas {
namespace grid {
namespace projection {

LonLatProjection::LonLatProjection(const eckit::Parametrisation& params) {
}

// copy constructor
LonLatProjection::LonLatProjection( const LonLatProjection& rhs ) : Projection(rhs) {
}

// clone method
Projection* LonLatProjection::clone() const  {
  return new LonLatProjection(*this);
}

eckit::geometry::LLPoint2 LonLatProjection::coords2lonlat(eckit::geometry::Point2 xy) const {

  return eckit::geometry::LLPoint2(xy[eckit::geometry::XX],xy[eckit::geometry::YY]);
}

eckit::geometry::Point2 LonLatProjection::lonlat2coords(eckit::geometry::LLPoint2 ll) const {

  return eckit::geometry::Point2(ll.lon(),ll.lat());
}

// specification
eckit::Properties LonLatProjection::spec() const {
  eckit::Properties proj_spec;
  proj_spec.set("projectionType",virtual_projection_type_str());
  return proj_spec;
}
register_BuilderT1(Projection,LonLatProjection,LonLatProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

