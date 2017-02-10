#include "atlas/grid/projection/RotatedSchmidtProjection.h"

#include "atlas/util/Constants.h"
#include <cmath>

#define D2R(X) (util::Constants::degreesToRadians()*(X))
#define R2D(X) (util::Constants::radiansToDegrees()*(X))

namespace atlas {
namespace grid {
namespace projection {

// constructor
RotatedSchmidtProjection::RotatedSchmidtProjection(const eckit::Parametrisation& params) :
  SchmidtProjection(params),
  Rotated(params) {
  // get stretching factor
  if( ! params.get("stretching_factor",c_) )
    throw eckit::BadParameter("stretching_factor missing in Params",Here());
}


// copy constructor
RotatedSchmidtProjection::RotatedSchmidtProjection( const RotatedSchmidtProjection& rhs ) :
    SchmidtProjection(rhs),
    Rotated(rhs) {
    c_=rhs.c_;
}

// clone method
Projection* RotatedSchmidtProjection::clone() const  {
  return new RotatedSchmidtProjection(*this);
}


eckit::geometry::LLPoint2 RotatedSchmidtProjection::coords2lonlat(eckit::geometry::Point2 xy) const {

  // stretch
  double lat=R2D(asin(cos(2*atan(1/c_*tan(acos(sin(D2R(xy[eckit::geometry::YY])))/2)))));

  eckit::geometry::LLPoint2 P(xy[eckit::geometry::XX],lat);

  // perform rotation
  rotate(P);

  return P;
}

eckit::geometry::Point2 RotatedSchmidtProjection::lonlat2coords(eckit::geometry::LLPoint2 P) const {

  // inverse rotation
  unrotate(P);

  // unstretch
  double lat2=R2D(asin(cos(2*atan(c_*tan(acos(sin(D2R(P.lat())))/2)))));

  return eckit::geometry::Point2(P.lon(),lat2);
}

// specification
eckit::Properties RotatedSchmidtProjection::spec() const {
  eckit::Properties proj_spec = Rotated::spec();
  proj_spec.set("projectionType",virtual_projection_type_str());
  proj_spec.set("projectionStretchingFactor",c_);
  return proj_spec;
}

register_BuilderT1(Projection,RotatedSchmidtProjection,RotatedSchmidtProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

