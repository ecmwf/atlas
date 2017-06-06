#include "atlas/domain/detail/Domain.h"
#include "atlas/projection/Projection.h"
#include "eckit/exception/Exceptions.h"

namespace atlas {
namespace domain {

Domain *Domain::create() {
  // default: global domain
  util::Config projParams;
  projParams.set("type","global");
  return Domain::create(projParams);
}

Domain *Domain::create(const eckit::Parametrisation &p) {

  std::string domain_type;
  if (p.get("type",domain_type)) {
    return eckit::Factory<Domain>::instance().get(domain_type).create(p);
  }

  // should return error here
  throw eckit::BadParameter("type missing in Params",Here());
}

bool Domain::includesNorthPole( const Projection& proj ) const {
  double north_pole[] = {0.,90.};
  proj.lonlat2xy(north_pole);
  return contains(north_pole[0],north_pole[1]);
}

bool Domain::includesSouthPole( const Projection& proj ) const {
  double south_pole[] = {0.,-90.};
  proj.lonlat2xy(south_pole);
  return contains(south_pole[0],south_pole[1]);
}

}  // namespace domain
}  // namespace atlas
