#include "atlas/grid/regular/RegularLonLat.h"

namespace atlas {
namespace grid {
namespace regular {

register_BuilderT1(Grid,RegularLonLat,RegularLonLat::grid_type_str());

std::string RegularLonLat::grid_type_str() {
    return "regular_lonlat";
}


std::string RegularLonLat::className() {
    return "atlas.grid.regular.RegularLonLat";
}

std::string RegularLonLat::shortName() const {
    std::ostringstream s;
    if ( nlonmin() == 2*nlat() && nlat()%2==0 ) {
      s << "L"<< nlat()/2;
    } else {
      s << "L"<< nlonmin() << "x" << nlat();
    }
    return s.str();
}

RegularLonLat::RegularLonLat(const util::Config& config) :
    GlobalLonLat()
{
    long nlon, nlat, N;

    // dimensions
    if ( config.get("N",N) ) {
      nlon=4*N;nlat=2*N+1;
    } else {
      if ( !config.get("nlon",nlon) || !config.get("nlat",nlat) ) {
        throw eckit::BadParameter("RegularLonLat requires either N, or (nlon,nlat)",Here());
      }
    }

    // default: no shift
    shiftLon_=false;
    shiftLat_=false;

    // perform setup
    setup(nlon, nlat);
}

RegularLonLat::RegularLonLat(long nlon, long nlat) : GlobalLonLat() {

    // default: no shift
    shiftLon_=false;
    shiftLat_=false;

    // perform setup
    GlobalLonLat::setup(nlon, nlat);
}

RegularLonLat::RegularLonLat(long N) : GlobalLonLat() {
    long nlon, nlat;

    // dimensions
    nlon=4*N;nlat=2*N+1;

    // default: no shift
    shiftLon_=false;
    shiftLat_=false;

    // perform setup
    GlobalLonLat::setup(nlon, nlat);
}

extern "C" {

    Structured* atlas__grid__regular__RegularLonLat(size_t nlon, size_t nlat) {
        return new RegularLonLat(nlon,nlat);
    }

}


}  // namespace regular
}  // namespace grid
}  // namespace atlas

