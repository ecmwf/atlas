#include "atlas/grid/domain/RectangularDomain.h"


namespace atlas {
namespace grid {
namespace domain {

RectangularDomain::RectangularDomain(const eckit::Parametrisation& params) {

	// check if bbox is present
	std::vector<double> v(4);
	if( ! params.get("bbox",v) )
    throw eckit::BadParameter("bbox missing in Params",Here());
   
  // store vector elements
	xmin_=v[0];
	xmax_=v[1];
	ymin_=v[2];
	ymax_=v[3];
	
	// setup
	setup();
}

void RectangularDomain::setup() {
	// normalize domain: make sure xmax>=xmin and ymax>=ymin
	double swp;

	if (xmin_>xmax_) {
		swp=xmin_;xmin_=xmax_;xmax_=swp;
	}

	if (ymin_>ymax_) {
		swp=ymin_;ymin_=ymax_;ymax_=swp;
	}
	
	//std::cout << domain_type_str() << std::endl;
}

bool RectangularDomain::contains(eckit::geometry::Point2 xy) {
	// probably should be done with some margin ...
	return ( xmin_ <= xy[eckit::geometry::XX] && xmax_ >= xy[eckit::geometry::XX] && ymin_ <= xy[eckit::geometry::YY] && ymax_ >= xy[eckit::geometry::YY] );
}

register_BuilderT1(Domain,RectangularDomain,RectangularDomain::domain_type_str());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

