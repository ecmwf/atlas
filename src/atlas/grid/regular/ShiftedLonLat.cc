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

std::string ShiftedLonLat::shortName() const {
    std::ostringstream s;
    long mlat=nlat()+1;	// for shifted lat, one latitude was removed.
    if ( nlonmin() == 2*mlat && mlat%2==0 ) {
    	s << "S"<< mlat/2;
    } else {
	    s << "S"<< nlonmin() << "x" << mlat;
	  }
    return s.str();
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

ShiftedLonLat::ShiftedLonLat(long nlon, long nlat) {
		// set shift
		shiftLon_=true;
		shiftLat_=true;
		
		// perform setup
		GlobalLonLat::setup(nlon,nlat);
}

}  // namespace regular
}  // namespace grid
}  // namespace atlas

