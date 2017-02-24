#include <cmath>
#include "atlas/grid/projection/Projection.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {
namespace projection {

Projection* Projection::create() {
  // default: no projection, i.e. stay in (lon,lat)-space
  util::Config projParams;
  projParams.set("type","lonlat");
  return Projection::create(projParams);
}

Projection* Projection::create(const eckit::Parametrisation& p) {
  std::string projectionType;
  if (p.get("type",projectionType)) {
    return eckit::Factory<Projection>::instance().get(projectionType).create(p);
  }

  // should return error here
  throw eckit::BadParameter("type missing in Params",Here());
}


}  // namespace projection
}  // namespace grid
}  // namespace atlas

