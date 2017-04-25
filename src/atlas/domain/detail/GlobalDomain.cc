#include "atlas/domain/detail/GlobalDomain.h"

namespace atlas {
namespace domain {

namespace {
  constexpr std::array<double,2> yrange() { 
    return { -90., 90. };
  }
}

GlobalDomain::GlobalDomain() :
  ZonalBandDomain( yrange() ) {
}

GlobalDomain::GlobalDomain(const eckit::Parametrisation& p) : 
  GlobalDomain() {
}

eckit::Properties GlobalDomain::spec() const {
  eckit::Properties domain_prop;
  domain_prop.set("type",type());
  return domain_prop;
}

void GlobalDomain::print(std::ostream& os) const {
  os << "GlobalDomain";
}

register_BuilderT1(Domain,GlobalDomain,GlobalDomain::static_type());

}  // namespace domain
}  // namespace atlas

