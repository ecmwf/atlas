#include "atlas/grid/domain/RectangularDomain.h"


namespace atlas {
namespace grid {
namespace domain {

namespace {
  
  static bool is_global(double xmin, double xmax, double ymin, double ymax) {
    const double eps = 1.e-12;
    return std::abs( (xmax-xmin) - 360. ) < eps 
        && std::abs( (ymax-ymin) - 180. ) < eps ;
  }
  
}

RectangularDomain::RectangularDomain(const eckit::Parametrisation& params) {

  if( ! params.get("xmin",xmin_) )
    throw eckit::BadParameter("xmin missing in Params",Here());

  if( ! params.get("xmax",xmax_) )
    throw eckit::BadParameter("xmax missing in Params",Here());

  if( ! params.get("ymin",ymin_) )
    throw eckit::BadParameter("ymin missing in Params",Here());

  if( ! params.get("ymax",ymax_) )
    throw eckit::BadParameter("ymax missing in Params",Here());

  // normalize domain: make sure xmax>=xmin and ymax>=ymin
  double swp;

  if (xmin_>xmax_) {
    swp=xmin_;xmin_=xmax_;xmax_=swp;
  }

  if (ymin_>ymax_) {
    swp=ymin_;ymin_=ymax_;ymax_=swp;
  }
  
  global_ = is_global(xmin_,xmax_,ymin_,ymax_);

}

RectangularDomain::RectangularDomain( const std::array<double,2>& xrange, const std::array<double,2>& yrange ) {
  xmin_ = xrange[0];
  xmax_ = xrange[1];
  ymin_ = yrange[0];
  ymax_ = yrange[1];
  global_ = is_global(xmin_,xmax_,ymin_,ymax_);
}


bool RectangularDomain::contains(double x, double y) const {
  // probably should be done with some margin ...
  return ( xmin_ <= x && xmax_ >= x && ymin_ <= y && ymax_ >= y );
}

eckit::Properties RectangularDomain::spec() const {
  eckit::Properties domain_prop;
  domain_prop.set("domainType",virtual_domain_type_str());
  domain_prop.set("xmin",xmin_);
  domain_prop.set("xmax",xmax_);
  domain_prop.set("ymin",ymin_);
  domain_prop.set("ymax",ymax_);
  return domain_prop;
}

register_BuilderT1(Domain,RectangularDomain,RectangularDomain::domain_type_str());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

