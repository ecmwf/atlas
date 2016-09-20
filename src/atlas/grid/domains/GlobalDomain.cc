#include "atlas/grid/domains/GlobalDomain.h"


namespace atlas {
namespace grid {

GlobalDomain::GlobalDomain() {
	setup();
}

void GlobalDomain::setup() {

}

bool GlobalDomain::contains(eckit::geometry::Point2 xy) {
	return true;
}

}  // namespace grid
}  // namespace atlas

