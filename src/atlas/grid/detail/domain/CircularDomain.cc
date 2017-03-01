#include "atlas/grid/detail/domain/CircularDomain.h"


namespace atlas {
namespace grid {
namespace domain {

CircularDomain::CircularDomain(const eckit::Parametrisation& params) {
  // read data from params
  std::vector<double> centre(2);
  if ( !params.get("centre",centre) ) throw eckit::BadParameter("center missing in Params",Here());
  xc_=centre[0];yc_=centre[1];

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
  domain_prop.set("type",type());
  domain_prop.set("radius",radius_);
  std::vector<double> centre(2);
  centre[0] = xc_;
  centre[1] = yc_;
  domain_prop.set("centre",eckit::makeVectorValue(centre));
  return domain_prop;
}

std::string CircularDomain::units() const {
  NOTIMP;
}

void CircularDomain::print(std::ostream& os) const {
  os << "CircularDomain[radius=" << radius_ << "," << "centre=(" << xc_ << "," << yc_ << ")]";
}


register_BuilderT1(Domain,CircularDomain,CircularDomain::static_type());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

