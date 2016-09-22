#include "atlas/grid/domain/GlobalDomain.h"


namespace atlas {
namespace grid {
namespace domain {

GlobalDomain::GlobalDomain(const eckit::Parametrisation& p) {
	setup();
}

void GlobalDomain::setup() {

}

bool GlobalDomain::contains(eckit::geometry::Point2 xy) {
	return true;
}

register_BuilderT1(Domain,GlobalDomain,GlobalDomain::domain_type_str());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

