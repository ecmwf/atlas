#include "eckit/utils/MD5.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/projection/detail/LonLatProjection.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace projection {
namespace detail {

ProjectionImpl* ProjectionImpl::create() {
  // default: no projection, i.e. stay in (lon,lat)-space
  return new LonLatProjection();
}

ProjectionImpl* ProjectionImpl::create(const eckit::Parametrisation& p) {
  std::string projectionType;
  if (p.get("type",projectionType)) {
    return eckit::Factory<ProjectionImpl>::instance().get(projectionType).create(p);
  }

  // should return error here
  throw eckit::BadParameter("type missing in Params",Here());
}


Rotated::Rotated( const PointLonLat& south_pole, double rotation_angle ) :
  util::Rotation(south_pole,rotation_angle) {
}

Rotated::Rotated(const eckit::Parametrisation& p) :
  util::Rotation(p){
}

void Rotated::spec(Spec& s) const {
  std::vector<double> npole{ northPole().lon(), northPole().lat() };
  std::vector<double> spole{ southPole().lon(), southPole().lat() };
  s.set("north_pole",npole);
  s.set("south_pole",spole);
  s.set("rotation_angle",rotationAngle());
}

void Rotated::hash( eckit::Hash& hsh ) const {
  hsh.add("rotated");
  hsh.add(southPole().lon());
  hsh.add(southPole().lat());
  hsh.add(rotationAngle());
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas

