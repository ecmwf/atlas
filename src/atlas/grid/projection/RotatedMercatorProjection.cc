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
RotatedMercatorProjection::RotatedMercatorProjection(const eckit::Parametrisation& params) :
  MercatorProjection(params),
  Rotated(params) {
  lon0_=0.0;
}

// copy constructor
RotatedMercatorProjection::RotatedMercatorProjection( const RotatedMercatorProjection& rhs ) :
  MercatorProjection(rhs),
  Rotated(rhs) {
  lon0_=rhs.lon0_;
}

// clone method
Projection* RotatedMercatorProjection::clone() const  {
  return new RotatedMercatorProjection(*this);
}

// projection
eckit::geometry::Point2 RotatedMercatorProjection::lonlat2coords(eckit::geometry::LLPoint2 ll) const {

  // first unrotate
  eckit::geometry::LLPoint2 P(ll.lon(),ll.lat());
  unrotate(P);

  // then project
  return MercatorProjection::lonlat2coords(P);
}

// inverse projection
eckit::geometry::LLPoint2 RotatedMercatorProjection::coords2lonlat(eckit::geometry::Point2 xy) const {

  // inverse projection
  eckit::geometry::LLPoint2 P=MercatorProjection::coords2lonlat(xy);

  // unrotate
  rotate(P);

  // then project
  return P;
}

// specification
eckit::Properties RotatedMercatorProjection::spec() const {
  eckit::Properties proj_spec = Rotated::spec();
  proj_spec.set("projectionType",virtual_projection_type_str());
  proj_spec.set("projectionRadius",radius_);
  return proj_spec;
}


register_BuilderT1(Projection,RotatedMercatorProjection,RotatedMercatorProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

