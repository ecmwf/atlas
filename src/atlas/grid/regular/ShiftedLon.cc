#include "atlas/grid/regular/ShiftedLon.h"


namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,ShiftedLon,ShiftedLon::grid_type_str());

std::string ShiftedLon::grid_type_str() {
    return "shiftedLon";
}

std::string ShiftedLon::className() {
    return "atlas.grid.regular.ShiftedLon";
}

std::string ShiftedLon::shortName() const {
    std::ostringstream s;
    if ( nlonmin() == 2*nlat() && nlat()%2==0 ) {
      s << "Slon"<< nlat()/2;
    } else {
      s << "Slon"<< nlonmin() << "x" << nlat();
    }
    return s.str();
}

ShiftedLon::ShiftedLon(const util::Config& config)
{
    long nlon, nlat, N;

    // dimensions
    if ( config.get("N",N) ) {
      nlon=4*N;nlat=2*N;
    } else {
      if ( !config.get("nlon",nlon) || !config.get("nlat",nlat) ) {
        throw eckit::BadParameter("ShiftedLon requires either N, or (nlon,nlat)",Here());
      }
    }

    // set shift
    shiftLon_=true;
    shiftLat_=false;

    // perform setup
    GlobalLonLat::setup(nlon,nlat);
}

ShiftedLon::ShiftedLon(long nlon, long nlat) {
    // set shift
    shiftLon_=true;
    shiftLat_=false;

    // perform setup
    GlobalLonLat::setup(nlon,nlat);
}

extern "C" {


    Structured* atlas__grid__regular__ShiftedLon(size_t nlon, size_t nlat) {
        return new ShiftedLon(nlon,nlat);
    }

}
}  // namespace regular
}  // namespace grid
}  // namespace atlas

