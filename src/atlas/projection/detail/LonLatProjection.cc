#include "eckit/utils/MD5.h"
#include "atlas/projection/detail/LonLatProjection.h"

namespace atlas {
namespace projection {
namespace detail {

template <typename Rotation>
LonLatProjectionT<Rotation>::LonLatProjectionT( const eckit::Parametrisation& config ) :
  ProjectionImpl(),
  rotation_(config) {
}

template <typename Rotation>
LonLatProjectionT<Rotation>::Spec LonLatProjectionT<Rotation>::spec() const {
  Spec proj_spec;
  proj_spec.set("type",static_type());
  rotation_.spec(proj_spec);
  return proj_spec;
}

template <typename Rotation>
void LonLatProjectionT<Rotation>::hash( eckit::Hash& hsh ) const {
  hsh.add(static_type());
  rotation_.hash(hsh);
}

register_BuilderT1(ProjectionImpl,LonLatProjection,LonLatProjection::static_type());
register_BuilderT1(ProjectionImpl,RotatedLonLatProjection,RotatedLonLatProjection::static_type());

}  // namespace detail
}  // namespace projection
}  // namespace atlas

