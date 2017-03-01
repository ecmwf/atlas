#include "atlas/grid/detail/grid/regular/ShiftedLonLat.h"


namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

register_BuilderT1(Grid,ShiftedLonLat,ShiftedLonLat::static_type());

std::string ShiftedLonLat::static_type() {
    return "shifted_lonlat";
}

std::string ShiftedLonLat::name() const {
    std::ostringstream s;
    size_t mlat=nlat()+1;  // for shifted lat, one latitude was removed.
    if ( nlonmin() == 2*mlat && mlat%2==0 ) {
      s << "S"<< mlat/2;
    } else {
      s << "S"<< nlonmin() << "x" << nlat();
    }
    return s.str();
}

ShiftedLonLat::ShiftedLonLat( const Config& config )
{
    size_t nlon, nlat, N;

    // dimensions
    if ( config.get("N",N) ) {
      nlon=4*N;nlat=2*N;
    } else {
      if ( !config.get("nx",nlon) || !config.get("ny",nlat) ) {
        throw eckit::BadParameter("ShiftedLonLat requires either N, or (nx,ny)",Here());
      }
    }

    // perform setup
    GlobalLonLat::setup(nlon,nlat,Shift(true,true));
}

ShiftedLonLat::ShiftedLonLat( long nlon, long nlat ) {

    // perform setup
    GlobalLonLat::setup(nlon,nlat,Shift(true,true));
}

extern "C" {


    Structured* atlas__grid__regular__ShiftedLonLat(size_t nlon, size_t nlat) {
        return new ShiftedLonLat(nlon,nlat);
    }
}

}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas

