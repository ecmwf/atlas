#include "atlas/grid/regular/ShiftedLat.h"


namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,ShiftedLat,ShiftedLat::grid_type_str());

std::string ShiftedLat::grid_type_str() {
    return "shiftedLat";
}


std::string ShiftedLat::className() {
    return "atlas.grid.regular.ShiftedLat";
}


ShiftedLat::ShiftedLat(const util::Config& config)
{
		long nlon, nlat, N;
		
		// dimensions
		if ( config.get("N",N) ) {
			nlon=4*N;nlat=2*N;
		} else {
			if ( !config.get("nlon",nlon) || !config.get("nlat",nlat) ) {
				throw eckit::BadParameter("ShiftedLat requires either N, or (nlon,nlat)",Here());
			}
		}
		
		// set shift
		shiftLat_=true;
		
		// perform setup
		GlobalLonLat::setup(nlon,nlat);
}

eckit::Properties ShiftedLat::spec() const {
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

