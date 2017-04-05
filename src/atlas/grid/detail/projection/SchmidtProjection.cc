#include <cmath>
#include "eckit/utils/MD5.h"
#include "atlas/grid/detail/projection/SchmidtProjection.h"
#include "atlas/util/Constants.h"

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

// constructor
template <typename Rotation>
SchmidtProjectionT<Rotation>::SchmidtProjectionT(const eckit::Parametrisation& params) :
  Projection(),
  rotation_(params) {
  if( ! params.get("stretching_factor",c_) )
    throw eckit::BadParameter("stretching_factor missing in Params",Here());
}

template <typename Rotation>
void SchmidtProjectionT<Rotation>::xy2lonlat(double crd[]) const {

  // stretch
  crd[1]=R2D(std::asin(std::cos(2.*std::atan(1/c_*std::tan(std::acos(std::sin(D2R(crd[1])))*0.5)))));

  // perform rotation
  rotation_.rotate(crd);
}

template <typename Rotation>
void SchmidtProjectionT<Rotation>::lonlat2xy(double crd[]) const {

  // inverse rotation
  rotation_.unrotate(crd);

  // unstretch
  crd[1]=R2D(std::asin(std::cos(2.*std::atan(c_*std::tan(std::acos(std::sin(D2R(crd[1])))*0.5)))));
}

// specification
template <typename Rotation>
eckit::Properties SchmidtProjectionT<Rotation>::spec() const {
  eckit::Properties proj_spec;
  proj_spec.set("type",static_type());
  proj_spec.set("stretching_factor",c_);
  rotation_.spec(proj_spec);
  return proj_spec;
}

template <typename Rotation>
void SchmidtProjectionT<Rotation>::hash( eckit::MD5& md5 ) const {
  md5.add(static_type());
  rotation_.hash(md5);
  md5.add(c_);
}

register_BuilderT1(Projection,SchmidtProjection,SchmidtProjection::static_type());
register_BuilderT1(Projection,RotatedSchmidtProjection,RotatedSchmidtProjection::static_type());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

