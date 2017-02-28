#include <utility>
#include "atlas/grid/detail/domain/RectangularDomain.h"


namespace atlas {
namespace grid {
namespace domain {

namespace {
  
  static bool is_global(double xmin, double xmax, double ymin, double ymax, const std::string& units) {
    return false;

    // if( units != "degrees" )
    //   return false;
    //
    // const double eps = 1.e-12;
    // return std::abs( (xmax-xmin) - 360. ) < eps
    //     && std::abs( (ymax-ymin) - 180. ) < eps ;
  }
  
  static std::array<double,2> get_interval_x(const eckit::Parametrisation& params) {
    double xmin, xmax;
  
    if( ! params.get("xmin",xmin) )
      throw eckit::BadParameter("xmin missing in Params",Here());

    if( ! params.get("xmax",xmax) )
      throw eckit::BadParameter("xmax missing in Params",Here());
  
    return {xmin,xmax};
  }

  static std::array<double,2> get_interval_y(const eckit::Parametrisation& params) {
    double ymin, ymax;
  
    if( ! params.get("ymin",ymin) )
      throw eckit::BadParameter("ymin missing in Params",Here());

    if( ! params.get("ymax",ymax) )
      throw eckit::BadParameter("ymax missing in Params",Here());
  
    return {ymin,ymax};
  }
  
  static std::string get_units(const eckit::Parametrisation& params) {
    std::string units;
    if( ! params.get("units",units) )
      throw eckit::BadParameter("units missing in Params",Here());
    return units;
  }

}

RectangularDomain::RectangularDomain(const eckit::Parametrisation& params) :
  RectangularDomain( get_interval_x(params), get_interval_y(params), get_units(params) ) {
}

RectangularDomain::RectangularDomain( const Interval& x, const Interval& y, const std::string& units ) :
  xmin_(x[0]),
  xmax_(x[1]),
  ymin_(y[0]),
  ymax_(y[1]),
  units_(units) {

  // Make sure xmax>=xmin and ymax>=ymin
  if (xmin_>xmax_) std::swap(xmin_,xmax_);
  if (ymin_>ymax_) std::swap(ymin_,ymax_);
  global_ = is_global(xmin_,xmax_,ymin_,ymax_,units_);
}


bool RectangularDomain::contains(double x, double y) const {
  // probably should be done with some margin ...
  return ( xmin_ <= x && xmax_ >= x && ymin_ <= y && ymax_ >= y );
}

eckit::Properties RectangularDomain::spec() const {
  eckit::Properties domain_prop;
  domain_prop.set("domainType",type());
  domain_prop.set("xmin",xmin_);
  domain_prop.set("xmax",xmax_);
  domain_prop.set("ymin",ymin_);
  domain_prop.set("ymax",ymax_);
  return domain_prop;
}

void RectangularDomain::print(std::ostream& os) const {
  os << "RectangularDomain["
     <<  "xmin=" << xmin()
     << ",xmax=" << xmax()
     << ",ymin=" << ymin()
     << ",ymax=" << ymax()
     << ",units=" << units()
     << "]";
}

register_BuilderT1(Domain,RectangularDomain,RectangularDomain::static_type());

}  // namespace domain
}  // namespace grid
}  // namespace atlas

