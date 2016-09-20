#include "atlas/grid/domains/EmptyDomain.h"


namespace atlas {
namespace grid {

EmptyDomain::EmptyDomain() {
	setup();
}

void EmptyDomain::setup() {

}

bool EmptyDomain::contains(eckit::geometry::Point2 xy) {
	return false;
}

}  // namespace grid
}  // namespace atlas

