#include "atlas/grid/reduced/ClassicGaussian.h"

#include "atlas/grid/spacing/Spacing.h"
#include "atlas/grid/reduced/pl/classic/PointsPerLatitude.h"

namespace atlas {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,ClassicGaussian,ClassicGaussian::grid_type_str());

std::string ClassicGaussian::grid_type_str() {
    return "classicGaussian";
}


std::string ClassicGaussian::className() {
    return "atlas.grid.reduced.ClassicGaussian";
}

std::string ClassicGaussian::shortName() const {
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
  projection_.reset( projection::Projection::create(config_proj) );

  // setup
  setup(N);
}

ClassicGaussian::ClassicGaussian(const util::Config& config) : ReducedGaussian() {
    size_t N;

    // get N from config
    if ( !config.get("N",N) ) throw eckit::BadParameter("ClassicGaussian requires N",Here());

    // projection is lonlat
    util::Config config_proj;
    config_proj.set("type","lonlat");
    projection_.reset( projection::Projection::create(config_proj) );

    // setup
    setup(N);
}


}  // namespace regular
}  // namespace grid
}  // namespace atlas

