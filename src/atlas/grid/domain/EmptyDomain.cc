#include "atlas/grid/domain/EmptyDomain.h"


namespace atlas {
namespace grid {
namespace domain {

EmptyDomain::EmptyDomain(const eckit::Parametrisation& p) {
	setup();
}

void EmptyDomain::setup() {

}

bool EmptyDomain::contains(eckit::geometry::Point2 xy) {
	return false;
}

register_BuilderT1(Domain,EmptyDomain,EmptyDomain::domain_type_str());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

