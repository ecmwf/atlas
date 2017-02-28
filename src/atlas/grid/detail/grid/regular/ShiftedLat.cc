#include "atlas/grid/detail/grid/regular/ShiftedLat.h"


namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

register_BuilderT1(Grid,ShiftedLat,ShiftedLat::grid_type_str());

std::string ShiftedLat::grid_type_str() {
    return "shifted_lat";
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

ShiftedLat::ShiftedLat( const Config& config )
{
    long nlon, nlat, N;

    // dimensions
    if ( config.get("N",N) ) {
      nlon=4*N;nlat=2*N;
    } else {
      if ( !config.get("nx",nlon) || !config.get("ny",nlat) ) {
        throw eckit::BadParameter("ShiftedLat requires either N, or (nx,ny)",Here());
      }
    }

    // perform setup
    GlobalLonLat::setup(nlon,nlat,Shift(false,true));
}

ShiftedLat::ShiftedLat(long nlon, long nlat) {

    // perform setup
    GlobalLonLat::setup(nlon,nlat,Shift(false,true));
}

extern "C" {


    Structured* atlas__grid__regular__ShiftedLat(size_t nlon, size_t nlat) {
        return new ShiftedLat(nlon,nlat);
    }
}


}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas

