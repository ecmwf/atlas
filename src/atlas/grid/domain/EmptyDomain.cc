#include "atlas/grid/domain/EmptyDomain.h"


namespace atlas {
namespace grid {
namespace domain {

EmptyDomain::EmptyDomain(const eckit::Parametrisation& p) {
}

eckit::Properties EmptyDomain::spec() const {
  eckit::Properties domain_prop;
  domain_prop.set("domainType",virtual_domain_type_str());
  return domain_prop;
}

register_BuilderT1(Domain,EmptyDomain,EmptyDomain::domain_type_str());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

