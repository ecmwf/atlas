#include "atlas/grid/domain/EmptyDomain.h"


namespace atlas {
namespace grid {
namespace domain {

EmptyDomain::EmptyDomain() {
}

EmptyDomain::EmptyDomain(const eckit::Parametrisation& p) {
}

eckit::Properties EmptyDomain::spec() const {
  eckit::Properties domain_prop;
  domain_prop.set("domainType",type());
  return domain_prop;
}

void EmptyDomain::print(std::ostream& os) const {
  os << "EmptyDomain";
}

std::string EmptyDomain::units() const {
  NOTIMP;
}

register_BuilderT1(Domain,EmptyDomain,EmptyDomain::static_type());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

