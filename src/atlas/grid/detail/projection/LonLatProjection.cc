#include "eckit/utils/MD5.h"
#include "atlas/grid/detail/projection/LonLatProjection.h"

namespace atlas {
namespace grid {
namespace projection {

template <typename Rotation>
LonLatProjectionT<Rotation>::LonLatProjectionT( const eckit::Parametrisation& config ) :
  Projection(),
  rotation_(config) {
}

template <typename Rotation>
eckit::Properties LonLatProjectionT<Rotation>::spec() const {
  eckit::Properties proj_spec;
  proj_spec.set("type",static_type());
  rotation_.spec(proj_spec);
  return proj_spec;
}

template <typename Rotation>
void LonLatProjectionT<Rotation>::hash( eckit::MD5& md5 ) const {
  md5.add(static_type());
  rotation_.hash(md5);
}


ShiftedLonLatProjection::ShiftedLonLatProjection( double lon, double lat ) :
  lon_(lon),
  lat_(lat) {
}
ShiftedLonLatProjection::ShiftedLonLatProjection( const eckit::Parametrisation& config ) {
  std::vector<double> shift{0.,0.};
  config.get("shift",shift);
  lon_ = shift[0];
  lat_ = shift[1];
}

eckit::Properties ShiftedLonLatProjection::spec() const {
  eckit::Properties proj_spec;
  proj_spec.set("type",static_type());
  std::vector<double> shift{lon_,lat_};
  proj_spec.set("shift",eckit::makeVectorValue(shift));
  return proj_spec;
}

void ShiftedLonLatProjection::hash( eckit::MD5& md5 ) const {
  md5.add(static_type());
  md5.add(lon_);
  md5.add(lat_);
}

register_BuilderT1(Projection,LonLatProjection,LonLatProjection::static_type());
register_BuilderT1(Projection,RotatedLonLatProjection,RotatedLonLatProjection::static_type());
register_BuilderT1(Projection,ShiftedLonLatProjection,ShiftedLonLatProjection::static_type());

}  // namespace projection
}  // namespace grid
}  // namespace atlas

