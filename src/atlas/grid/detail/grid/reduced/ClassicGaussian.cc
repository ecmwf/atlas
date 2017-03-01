#include "atlas/grid/detail/grid/reduced/ClassicGaussian.h"

#include "atlas/grid/detail/spacing/Spacing.h"
#include "atlas/grid/detail/grid/reduced/pl/classic/PointsPerLatitude.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,ClassicGaussian,ClassicGaussian::static_type());

std::string ClassicGaussian::static_type() {
    return "classic_gaussian";
}

std::string ClassicGaussian::name() const {
   std::ostringstream s;
   s << "N" << nlat()/2;
   return s.str();
}


void ClassicGaussian::setup(size_t N) {

    util::Config config_spacing, config_domain, config_proj;

    // number of longitudes: from predefined sets
    std::vector<long> pl(N);      // number of longitudes per latitude
    pl::classic::points_per_latitude_npole_equator(N,pl.data());

    // setup from reducedGaussian
    ReducedGaussian::setup(N,pl.data());

}

ClassicGaussian::ClassicGaussian(size_t N) : ReducedGaussian() {

  // projection is lonlat
  util::Config config_proj;
  config_proj.set("type","lonlat");
  projection_ = Projection(config_proj);

  // setup
  setup(N);
}

ClassicGaussian::ClassicGaussian(const util::Config& config) : ReducedGaussian() {
    size_t N;

    // get N from config
    if ( !config.get("N",N) ) throw eckit::BadParameter("ClassicGaussian requires N",Here());

    // default projection is lonlat
    util::Config config_proj;
    if( not config.get("projection",config_proj) ) {
      config_proj.set("type","lonlat");
    }
    projection_ = Projection(config_proj);

    // setup
    setup(N);
}


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas

