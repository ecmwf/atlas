#include "atlas/grid/domain/CircularDomain.h"


namespace atlas {
namespace grid {
namespace domain {

CircularDomain::CircularDomain(const eckit::Parametrisation& params) {
  // read data from params
  std::vector<double> xyc(2);
  if ( !params.get("center",xyc) ) throw eckit::BadParameter("center missing in Params",Here());

  xc_=xyc[0];yc_=xyc[1];

  if ( !params.get("radius",radius_) ) throw eckit::BadParameter("radius missing in Params",Here());

  setup();
}

void CircularDomain::setup() {

}

bool CircularDomain::contains(eckit::geometry::Point2 xy) const {
  // probably should be done with some margin ...
  return ( (xy[eckit::geometry::XX]-xc_)*(xy[eckit::geometry::XX]-xc_)+(xy[eckit::geometry::YY]-yc_)*(xy[eckit::geometry::YY]-yc_) <= radius_*radius_ );
}

eckit::Properties CircularDomain::spec() const {
  eckit::Properties domain_prop;
  domain_prop.set("domainType",virtual_domain_type_str());
  domain_prop.set("radius",radius_);
  std::vector<double> xyc(2);
  xyc[0] = xc_;
  xyc[1] = yc_;
  domain_prop.set("center",eckit::makeVectorValue(xyc));
  return domain_prop;
}

register_BuilderT1(Domain,CircularDomain,CircularDomain::domain_type_str());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

