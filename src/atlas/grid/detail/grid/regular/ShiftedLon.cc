#include "atlas/grid/detail/grid/regular/ShiftedLon.h"


namespace atlas {
namespace grid {
namespace detail {
namespace grid {
namespace regular {

register_BuilderT1(Grid,ShiftedLon,ShiftedLon::static_type());

std::string ShiftedLon::static_type() {
    return "shifted_lon";
}

std::string ShiftedLon::name() const {
    std::ostringstream s;
    if ( nlonmin() == 2*nlat() && nlat()%2==0 ) {
      s << "Slon"<< nlat()/2;
    } else {
      s << "Slon"<< nlonmin() << "x" << nlat();
    }
    return s.str();
}

ShiftedLon::ShiftedLon( const Config& config )
{
    long nlon, nlat, N;

    // dimensions
    if ( config.get("N",N) ) {
      nlon=4*N;nlat=2*N;
    } else {
      if ( !config.get("nx",nlon) || !config.get("ny",nlat) ) {
        throw eckit::BadParameter("ShiftedLon requires either N, or (nx,ny)",Here());
      }
    }

    // perform setup
    GlobalLonLat::setup(nlon,nlat,Shift(true,false));
}

ShiftedLon::ShiftedLon(long nlon, long nlat) {
    // perform setup
    GlobalLonLat::setup(nlon,nlat,Shift(true,false));
}

extern "C" {


    Structured* atlas__grid__regular__ShiftedLon(size_t nlon, size_t nlat) {
        return new ShiftedLon(nlon,nlat);
    }

}
}  // namespace regular
}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas

