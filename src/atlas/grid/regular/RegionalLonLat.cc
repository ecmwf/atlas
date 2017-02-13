#include "atlas/grid/regular/RegionalLonLat.h"


namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,RegionalLonLat,RegionalLonLat::grid_type_str());

std::string RegionalLonLat::grid_type_str() {
    return "regional_lonlat";
}


std::string RegionalLonLat::className() {
    return "atlas.grid.regular.RegionalLonLat";
}

void RegionalLonLat::setup(const util::Config& config) {

    util::Config config_proj, config_dom, config_rr;
    long nlon, nlat;
    std::vector<double> bbox(4);
    double east, west, south, north;

    // get domain boundaries
    if ( config.get("bounding_box",bbox) ) {
      east=bbox[0];
      west=bbox[1];
      south=bbox[2];
      north=bbox[3];
    } else {
      if ( ! ( config.get("east",east) && config.get("west",west) && config.get("south",south) && config.get("north",north) ) ) {
        throw eckit::BadParameter("RegionalLonLat grid domain should be specified by bbox, or by east, west, south and north.",Here());
      }
    }

    // perform checks on bounds
    if (south>north)
      throw eckit::BadValue("RegionalLonLat grid domain requires north>south.",Here());
    while (east<west) east+=360.;
    bbox[0]=east;
    bbox[1]=west;
    bbox[2]=south;
    bbox[3]=north;

    // define domain subconfiguration
    config_dom.set("bounding_box",bbox);
    config_dom.set("type","rectangular");

    // grid dimensions
    if ( !config.get("nlon",nlon) )
      throw eckit::BadParameter("RegionalLonLat grid domain requires nlon.",Here());
    if ( !config.get("nlat",nlat) )
      throw eckit::BadParameter("RegionalLonLat grid domain requires nlat.",Here());

    // projection is lonlat
    config_proj.set("type","lonlat");

    // put subconfigurations in master grid configuration for RegularRegional
    config_rr.set("domain",config_dom);
    config_rr.set("projection",config_proj);
    config_rr.set("nx",nlon);
    config_rr.set("ny",nlat);

    // setup a regular regional grid
    RegularRegional::setup(config_rr);

}

RegionalLonLat::RegionalLonLat(const util::Config& config) :
    RegularRegional()
{
    // perform setup
    setup(config);
}

RegionalLonLat::RegionalLonLat() : RegularRegional() {
}

std::string RegionalLonLat::shortName() const {
    std::ostringstream s;
    s << "RLL"<< nlonmin() << "x" << nlat();
    return s.str();
}

}  // namespace regular
}  // namespace grid
}  // namespace atlas

