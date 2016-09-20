#include "atlas/grid/domains/CircularDomain.h"


namespace atlas {
namespace grid {

CircularDomain::CircularDomain(double xc, double yc, double radius):
		xc_(xc), yc_(yc), radius_(radius) {
	setup();
}

	
CircularDomain::CircularDomain(eckit::geometry::Point2 XC, double radius ): radius_(radius) {
	xc_=XC[eckit::geometry::XX];
	yc_=XC[eckit::geometry::YY];
	
	setup();
}

void CircularDomain::setup() {

}

bool CircularDomain::contains(eckit::geometry::Point2 xy) {
	// probably should be done with some margin ...
	return ( (xy[eckit::geometry::XX]-xc_)*(xy[eckit::geometry::XX]-xc_)+(xy[eckit::geometry::YY]-yc_)*(xy[eckit::geometry::YY]-yc_) <= radius_*radius_ );
}

}  // namespace grid
}  // namespace atlas

