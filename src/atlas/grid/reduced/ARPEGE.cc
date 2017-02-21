#include "atlas/grid/reduced/ARPEGE.h"

#include "atlas/grid/projection/Projection.h"

namespace atlas {
namespace grid {
namespace reduced {

register_BuilderT1(Grid,ARPEGE,ARPEGE::grid_type_str());

std::string ARPEGE::grid_type_str() {
    return "ARPEGE";
}


std::string ARPEGE::className() {
    return "atlas.grid.reduced.ARPEGE";
}

std::string ARPEGE::shortName() const {
   std::ostringstream s;
   s << "ARP" << nlat()/2;
   return s.str();
}

ARPEGE::ARPEGE(const util::Config& config) :
    ClassicGaussian()
{
    size_t N;
    double c;
    util::Config config_proj;

    // get N from config
    if ( !config.get("N",N) ) throw eckit::BadParameter("ARPEGE requires N",Here());

    // set projection
    // get stretching factor; defaults to 1
    if ( !config.get("stretching_factor",c) ) c=1.0;
    config_proj.set("stretching_factor",c);
    std::vector<double> pole(2);
    if ( config.get("north_pole",pole) ) {
      config_proj.set("type","rotated_schmidt");
      config_proj.set("north_pole",pole);
    } else {
      config_proj.set("type","schmidt");
    }
    projection_.reset( projection::Projection::create(config_proj) );

    // setup
    ClassicGaussian::setup(N);
}

}  // namespace regular
}  // namespace grid
}  // namespace atlas

