
#include "atlas/grid/domain/Domain.h"
#include "atlas/util/Config.h"
#include "eckit/exception/Exceptions.h"

#warning TODO: Non-global grids may also include north-pole, south-pole, or be periodic-east-west

namespace atlas {
namespace grid {
namespace domain {

Domain *Domain::create() {
  // default: global domain
  util::Config projParams;
  projParams.set("domainType","global");
  return Domain::create(projParams);
}

Domain *Domain::create(const eckit::Parametrisation &p) {

  std::string domainType;
  if (p.get("domainType",domainType)) {
    return eckit::Factory<Domain>::instance().get(domainType).create(p);
  }

  // should return error here
  throw eckit::BadParameter("domainType missing in Params",Here());
}


}  // namespace domain
}  // namespace grid
}  // namespace atlas
