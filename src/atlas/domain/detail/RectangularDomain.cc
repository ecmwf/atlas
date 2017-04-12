#include <utility>
#include "atlas/domain/detail/RectangularDomain.h"


namespace atlas {
namespace domain {

using Interval = RectangularDomain::Interval;

namespace {

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

bool RectangularDomain::is_global(
   const Interval& x,
   const Interval& y,
   const std::string& units ) {

  if( units != "degrees" )
    return false;

  const double eps = 1.e-12;
  return std::abs( (x[1]-x[0]) - 360. ) < eps
      && std::abs( (y[1]-y[0]) - 180. ) < eps ;
}

bool RectangularDomain::is_zonal_band(
  const Interval& x,
  const std::string& units) {

  if( units != "degrees" )
    return false;

  const double eps = 1.e-12;
  return std::abs( (x[1]-x[0]) - 360. ) < eps;
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
  global_ = is_global( {xmin_,xmax_}, {ymin_,ymax_} ,units_);
  
  const double tol = 1.e-6;
  xmin_tol_ = xmin_-tol;
  ymin_tol_ = ymin_-tol;
  xmax_tol_ = xmax_+tol;
  ymax_tol_ = ymax_+tol;
}


bool RectangularDomain::contains(double x, double y) const {
  return contains_x(x) and contains_y(y);
}

eckit::Properties RectangularDomain::spec() const {
  eckit::Properties domain_prop;
  domain_prop.set("type",type());
  domain_prop.set("xmin",xmin());
  domain_prop.set("xmax",xmax());
  domain_prop.set("ymin",ymin());
  domain_prop.set("ymax",ymax());
  domain_prop.set("units",units());
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
}  // namespace atlas

