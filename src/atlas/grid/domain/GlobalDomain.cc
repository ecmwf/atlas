#include "atlas/grid/domain/GlobalDomain.h"


namespace atlas {
namespace grid {
namespace domain {

GlobalDomain::GlobalDomain(const eckit::Parametrisation& p) {
  setup();
}

void GlobalDomain::setup() {

}

bool GlobalDomain::contains(eckit::geometry::Point2 xy) const {
  return true;
}

eckit::Properties GlobalDomain::spec() const {
  eckit::Properties domain_prop;
  domain_prop.set("domainType",virtual_domain_type_str());
  return domain_prop;
}

register_BuilderT1(Domain,GlobalDomain,GlobalDomain::domain_type_str());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

