#include <cmath>
#include "eckit/utils/MD5.h"
#include "atlas/projection/detail/MercatorProjection.h"
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
namespace projection {
namespace detail {

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

template <typename Rotation>
void MercatorProjectionT<Rotation>::lonlat2xy(double crd[]) const {

  // first unrotate
  rotation_.unrotate(crd);

  // then project
  crd[0] = radius_*(D2R(crd[0]-lon0_));
  crd[1] = radius_*std::log(std::tan(D2R(45.+crd[1]*0.5)));
}

template <typename Rotation>
void MercatorProjectionT<Rotation>::xy2lonlat(double crd[]) const {

  // first projection
  crd[0] = lon0_ + R2D(crd[0]*inv_radius_);
  crd[1] = 2.*R2D(std::atan(std::exp(crd[1]*inv_radius_)))-90.;

  // then rotate
  rotation_.rotate(crd);
}


// specification
template <typename Rotation>
eckit::Properties MercatorProjectionT<Rotation>::spec() const {
  eckit::Properties proj_spec;
  proj_spec.set("type",static_type());
  proj_spec.set("longitude0",lon0_);
  if( radius_ != util::Earth::radiusInMeters() ) proj_spec.set("radius",radius_);
  rotation_.spec(proj_spec);
  return proj_spec;
}

template <typename Rotation>
void MercatorProjectionT<Rotation>::hash( eckit::MD5& md5 ) const {
  md5.add(static_type());
  rotation_.hash(md5);
  md5.add(lon0_);
  md5.add(radius_);
}

register_BuilderT1(Projection,MercatorProjection,MercatorProjection::static_type());
register_BuilderT1(Projection,RotatedMercatorProjection,RotatedMercatorProjection::static_type());

}  // namespace detail
}  // namespace projection
}  // namespace atlas

