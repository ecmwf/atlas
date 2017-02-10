#include <cmath>
#include "atlas/grid/projection/MercatorProjection.h"
#include "atlas/util/Constants.h"

/*
Projection formula's for Mercator projection from "Map Projections: A Working Manual"

The origin of the xy-system is at (lon0,0)

*/

namespace {
  static double D2R(const double x) {
    return atlas::util::Constants::degreesToRadians()*x;
  }
  static double R2D(const double x) {
    return atlas::util::Constants::radiansToDegrees()*x;
  }
}

namespace atlas {
namespace grid {
namespace projection {

// constructors
template <typename Rotation>
MercatorProjectionT<Rotation>::MercatorProjectionT(const eckit::Parametrisation& params) :
  Projection(),
  rotation_(params) {

  // check presence of radius
  if( ! params.get("radius",radius_) )
    radius_=util::Earth::radiusInMeters();
  // check presence of lon0
  if( ! params.get("longitude0",lon0_) )
    lon0_=0.0;

}

// copy constructor
template <typename Rotation>
MercatorProjectionT<Rotation>::MercatorProjectionT( const MercatorProjectionT& rhs ) :
  Projection(rhs),
  rotation_(rhs.rotation_) {
  radius_=rhs.radius_;
  lon0_=rhs.lon0_;
}

// clone method
template <typename Rotation>
Projection* MercatorProjectionT<Rotation>::clone() const  {
  return new MercatorProjectionT<Rotation>(*this);
}

// projection
template <typename Rotation>
eckit::geometry::Point2 MercatorProjectionT<Rotation>::lonlat2coords(eckit::geometry::LLPoint2 ll) const {
  // first unrotate
  rotation_.unrotate(ll);

  // then project
  double x=radius_*(D2R(ll.lon()-lon0_));
  double y=radius_*std::log(std::tan(D2R(45.+ll.lat()*0.5)));

  return eckit::geometry::Point2(x,y);
}

// inverse projection
template <typename Rotation>
eckit::geometry::LLPoint2 MercatorProjectionT<Rotation>::coords2lonlat(eckit::geometry::Point2 coords) const {

  // first projection
  const double x=coords[eckit::geometry::XX];
  const double y=coords[eckit::geometry::YY];
  eckit::geometry::LLPoint2 lonlat(
      lon0_+R2D(x/radius_),
      2.*R2D(std::atan(std::exp(y/radius_)))-90.
  );

  // then rotate
  rotation_.rotate(lonlat);

  // then project
  return lonlat;
}

// specification
template <typename Rotation>
eckit::Properties MercatorProjectionT<Rotation>::spec() const {
  eckit::Properties proj_spec;
  proj_spec.set("projectionType",virtual_projection_type_str());
  proj_spec.set("projectionLongitude0",lon0_);
  proj_spec.set("projectionRadius",radius_);
  rotation_.spec(proj_spec);
  return proj_spec;
}

register_BuilderT1(Projection,MercatorProjection,MercatorProjection::projection_type_str());
register_BuilderT1(Projection,RotatedMercatorProjection,RotatedMercatorProjection::projection_type_str());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

