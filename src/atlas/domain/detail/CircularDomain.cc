#include "atlas/domain/detail/CircularDomain.h"


namespace atlas {
namespace domain {

CircularDomain::CircularDomain(const std::array<double,2>& centre, const double& radius, const std::string& units ) :
  xc_( centre[0] ),
  yc_( centre[1] ),
  radius_( radius ),
  rr_( radius*radius ),
  units_(units) {
}


CircularDomain::CircularDomain(const eckit::Parametrisation& params) {
  // read data from params
  std::vector<double> centre(2);
  if ( !params.get("centre",centre) ) throw eckit::BadParameter("center missing in Params",Here());
  xc_=centre[0];yc_=centre[1];

  if ( !params.get("radius",radius_) ) throw eckit::BadParameter("radius missing in Params",Here());

  if ( !params.get("units",units_) ) throw eckit::BadParameter("units missing in Params",Here());

  rr_ = radius_*radius_;
}


bool CircularDomain::contains(double x, double y) const {
  // probably should be done with some margin ...
  double xx = x-xc_; xx *= xx;
  double yy = y-yc_; yy *= yy;
  return ( xx+yy <= rr_ );
}

CircularDomain::Spec CircularDomain::spec() const {
  Spec domain_spec;
  domain_spec.set("type",type());
  domain_spec.set("radius",radius_);
  domain_spec.set("centre", std::vector<double>{xc_,yc_} );
  return domain_spec;
}

void CircularDomain::hash(eckit::Hash& h) const {
  spec().hash(h);
}

std::string CircularDomain::units() const {
  return units_;
}

void CircularDomain::print(std::ostream& os) const {
  os << "CircularDomain[radius=" << radius_ << "," << "centre=(" << xc_ << "," << yc_ << ")]";
}


register_BuilderT1(Domain,CircularDomain,CircularDomain::static_type());

}  // namespace domain
}  // namespace atlas

