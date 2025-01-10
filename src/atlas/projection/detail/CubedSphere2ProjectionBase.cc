#include "CubedSphere2ProjectionBase.h"

#include "atlas/runtime/Trace.h"
#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"
#include "eckit/utils/Hash.h"

namespace atlas {
namespace projection {
namespace detail {

using Spec = atlas::util::Config;

CubedSphere2ProjectionBase::CubedSphere2ProjectionBase(const eckit::Parametrisation& params) {
    ATLAS_TRACE("CubedSphere2ProjectionBase::CubedSphere2ProjectionBase");
}

Spec CubedSphere2ProjectionBase::spec() const {
    Spec proj;
    proj.set("type", static_type());
    return proj;
}

Jacobian CubedSphere2ProjectionBase::jacobian(const PointLonLat& lonlat) const {
    ATLAS_NOTIMPLEMENTED;
}

void CubedSphere2ProjectionBase::hash(eckit::Hash& h) const {
    // Add to hash
    h.add(static_type());
}

namespace {
static ProjectionBuilder<CubedSphere2ProjectionBase> register_1(CubedSphere2ProjectionBase::static_type());
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
