#include "eckit/utils/MD5.h"
#include "atlas/grid/detail/projection/Projection.h"
#include "atlas/grid/detail/projection/LonLatProjection.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {
namespace projection {

Projection* Projection::create() {
  // default: no projection, i.e. stay in (lon,lat)-space
  return new LonLatProjection();
}

Projection* Projection::create(const eckit::Parametrisation& p) {
  std::string projectionType;
  if (p.get("type",projectionType)) {
    return eckit::Factory<Projection>::instance().get(projectionType).create(p);
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

void Rotated::spec(eckit::Properties& s) const {
  std::vector<double> npole{ northPole().lon(), northPole().lat() };
  std::vector<double> spole{ southPole().lon(), southPole().lat() };
  s.set("north_pole",eckit::makeVectorValue(npole));
  s.set("south_pole",eckit::makeVectorValue(spole));
  s.set("rotation_angle",rotationAngle());
}

void Rotated::hash( eckit::MD5& md5 ) const {
  md5.add("rotated");
  md5.add(southPole().lon());
  md5.add(southPole().lat());
  md5.add(rotationAngle());
}

}  // namespace projection
}  // namespace grid
}  // namespace atlas

