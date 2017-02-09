#include "atlas/grid/domain/CircularDomain.h"


namespace atlas {
namespace grid {
namespace domain {

CircularDomain::CircularDomain(const eckit::Parametrisation& params) {
  // read data from params
  std::vector<double> center(2);
  if ( !params.get("center",center) ) throw eckit::BadParameter("center missing in Params",Here());
  xc_=center[0];yc_=center[1];

  if ( !params.get("radius",radius_) ) throw eckit::BadParameter("radius missing in Params",Here());

  rr_ = radius_*radius_;
}


bool CircularDomain::contains(double x, double y) const {
  // probably should be done with some margin ...
  double xx = x-xc_; xx *= xx;
  double yy = y-yc_; yy *= yy;
  return ( xx+yy <= rr_ );
}

eckit::Properties CircularDomain::spec() const {
  eckit::Properties domain_prop;
  domain_prop.set("domainType",virtual_domain_type_str());
  domain_prop.set("radius",radius_);
  std::vector<double> center(2);
  center[0] = xc_;
  center[1] = yc_;
  domain_prop.set("center",eckit::makeVectorValue(center));
  return domain_prop;
}

register_BuilderT1(Domain,CircularDomain,CircularDomain::domain_type_str());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

