#include "atlas/grid/detail/grid/regular/ShiftedLat.h"


namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

register_BuilderT1(Grid,ShiftedLat,ShiftedLat::static_type());

std::string ShiftedLat::static_type() {
    return "shifted_lat";
}

std::string ShiftedLat::name() const {
    std::ostringstream s;
    size_t mlat=nlat()+1;  // for shifted lat, one latitude was removed.
    if ( nlonmin() == 2*mlat && mlat%2==0 ) {
      s << "Slat"<< mlat/2;
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

