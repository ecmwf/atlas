#include "atlas/grid/domain/RectangularDomain.h"


namespace atlas {
namespace grid {
namespace domain {

RectangularDomain::RectangularDomain(const eckit::Parametrisation& params) {

  // check if bbox is present
  std::vector<double> v(4);
  if( ! params.get("bounding_box",v) )
    throw eckit::BadParameter("bounding_box missing in Params",Here());

  // store vector elements
  xmin_=v[0];
  xmax_=v[1];
  ymin_=v[2];
  ymax_=v[3];

  // normalize domain: make sure xmax>=xmin and ymax>=ymin
  double swp;

  if (xmin_>xmax_) {
    swp=xmin_;xmin_=xmax_;xmax_=swp;
  }

  if (ymin_>ymax_) {
    swp=ymin_;ymin_=ymax_;ymax_=swp;
  }

  periodic_x_ = false;
  periodic_y_ = false;
  params.get("periodic_x",periodic_x_);
  params.get("periodic_y",periodic_y_);
}

bool RectangularDomain::contains(double x, double y) const {
  // probably should be done with some margin ...
  return ( xmin_ <= x && xmax_ >= x && ymin_ <= y && ymax_ >= y );
}

std::vector<double> RectangularDomain::bbox() const {

  std::vector<double> bb(4);
  bb[0]=xmin_;
  bb[1]=xmax_;
  bb[2]=ymin_;
  bb[3]=ymax_;

  return bb;

}

eckit::Properties RectangularDomain::spec() const {
  eckit::Properties domain_prop;
  domain_prop.set("domainType",virtual_domain_type_str());
  std::vector<double> bb(4);
  bb=bbox();
  domain_prop.set("bounding_box",eckit::makeVectorValue(bb));
  return domain_prop;
}

register_BuilderT1(Domain,RectangularDomain,RectangularDomain::domain_type_str());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

