#include "atlas/grid/projection/RotatedLonLatProjection.h"

#include "atlas/util/Constants.h"
#include <cmath>

#define D2R(X) (util::Constants::degreesToRadians()*(X))
#define R2D(X) (util::Constants::radiansToDegrees()*(X))

namespace atlas {
namespace grid {
namespace projection {

// constructor
RotatedLonLatProjection::RotatedLonLatProjection(const eckit::Parametrisation& params) :
  LonLatProjection(params),
  Rotated(params) {
}

// copy constructor
RotatedLonLatProjection::RotatedLonLatProjection( const RotatedLonLatProjection& rhs ) :
  LonLatProjection(rhs),
  Rotated(rhs) {
}

// clone method
Projection* RotatedLonLatProjection::clone() const  {
  return new RotatedLonLatProjection(*this);
}


eckit::geometry::LLPoint2 RotatedLonLatProjection::coords2lonlat(eckit::geometry::Point2 xy) const {

  // point
  eckit::geometry::LLPoint2 P(xy[eckit::geometry::XX],xy[eckit::geometry::YY]);

  // perform rotation
  rotate(P);

  return P;
}

eckit::geometry::Point2 RotatedLonLatProjection::lonlat2coords(eckit::geometry::LLPoint2 P) const {

  // inverse rotation
  unrotate(P);

  return eckit::geometry::Point2(P.lon(),P.lat());
}


// specification
eckit::Properties RotatedLonLatProjection::spec() const {
  eckit::Properties proj_spec = Rotated::spec();
  proj_spec.set("projectionType",virtual_projection_type_str());
  return proj_spec;
}

register_BuilderT1(Projection,RotatedLonLatProjection,RotatedLonLatProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

