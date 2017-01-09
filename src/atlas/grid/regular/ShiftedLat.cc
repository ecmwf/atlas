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

std::string ShiftedLat::shortName() const {
    std::ostringstream s;
    if ( nlonmin() == 2*nlat() ) {
    	s << "Slat"<< nlat()/2;
    } else {
	    s << "Slat"<< nlonmin() << "x" << nlat();
	  }
    return s.str();
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
		shiftLon_=false;
		shiftLat_=true;
		
		// perform setup
		GlobalLonLat::setup(nlon,nlat);
}

ShiftedLat::ShiftedLat(long nlon, long nlat) {
		// set shift
		shiftLon_=false;
		shiftLat_=true;
		
		// perform setup
		GlobalLonLat::setup(nlon,nlat);
}

}  // namespace regular
}  // namespace grid
}  // namespace atlas

