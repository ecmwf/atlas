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

  inv_radius_ = 1./radius_;
}

// copy constructor
template <typename Rotation>
MercatorProjectionT<Rotation>::MercatorProjectionT( const MercatorProjectionT& rhs ) :
  Projection(rhs),
  rotation_(rhs.rotation_) {
  radius_=rhs.radius_;
  lon0_=rhs.lon0_;
  inv_radius_ = 1./radius_;
}

// clone method
template <typename Rotation>
Projection* MercatorProjectionT<Rotation>::clone() const  {
  return new MercatorProjectionT<Rotation>(*this);
}


template <typename Rotation>
void MercatorProjectionT<Rotation>::lonlat2coords(double crd[]) const {

  // first unrotate
  rotation_.unrotate(crd);

  // then project
  crd[0] = radius_*(D2R(crd[0]-lon0_));
  crd[1] = radius_*std::log(std::tan(D2R(45.+crd[1]*0.5)));
}

template <typename Rotation>
void MercatorProjectionT<Rotation>::coords2lonlat(double crd[]) const {

  // first projection
  crd[0] = lon0_ + R2D(crd[0]*inv_radius_);
  crd[1] = 2.*R2D(std::atan(std::exp(crd[0]*inv_radius_)))-90.;

  // then rotate
  rotation_.rotate(crd);
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

