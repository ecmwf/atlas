#include "eckit/utils/MD5.h"
#include "atlas/projection/detail/LonLatProjection.h"

namespace atlas {
namespace projection {
namespace detail {

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

register_BuilderT1(Projection,LonLatProjection,LonLatProjection::static_type());
register_BuilderT1(Projection,RotatedLonLatProjection,RotatedLonLatProjection::static_type());

}  // namespace detail
}  // namespace projection
}  // namespace atlas

