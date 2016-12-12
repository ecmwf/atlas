#include "atlas/grid/regular/ShiftedLonLat.h"


namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,ShiftedLonLat,ShiftedLonLat::grid_type_str());

std::string ShiftedLonLat::grid_type_str() {
    return "shiftedLonLat";
}


std::string ShiftedLonLat::className() {
    return "atlas.grid.regular.ShiftedLonLat";
}


ShiftedLonLat::ShiftedLonLat(const util::Config& config)
{
		long nlon, nlat, N;
		
		// dimensions
		if ( config.get("N",N) ) {
			nlon=4*N;nlat=2*N;
		} else {
			if ( !config.get("nlon",nlon) || !config.get("nlat",nlat) ) {
				throw eckit::BadParameter("ShiftedLonLat requires either N, or (nlon,nlat)",Here());
			}
		}
		
		// set shift
		shiftLon_=true;
		shiftLat_=true;
		
		// perform setup
		GlobalLonLat::setup(nlon,nlat);
}

eckit::Properties ShiftedLonLat::spec() const {
    eckit::Properties grid_spec;
    grid_spec.set("grid_type",  gridType());
    /*
    grid_spec.set("short_name", shortName());
    grid_spec.set("N",    N());
    grid_spec.set("nlon", nlon());
    grid_spec.set("nlat", nlat());
    grid_spec.set("domain", domain_spec(domain_) );
    */
    return grid_spec;
}


}  // namespace regular
}  // namespace grid
}  // namespace atlas

