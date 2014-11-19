// TL47

#include "atlas/grids/reduced_gg/Grids.h"

namespace atlas {
namespace grids {
namespace reduced_gg {

N24::regist N24_builders;

void N24::construct()
{
  int N=24;
  int lon[] = {
    20,
    25,
    36,
    40,
    45,
    48,
    54,
    60,
    64,
    72,
    80,
    80,
    90,
    90,
    96,
    96,
    96,
    96,
    96,
    96,
    96,
    96,
    96,
    96
  };
  std::vector<double> lat(N);
  eckit::Log::warning() << className() << " uses predicted gaussian latitudes" << std::endl;
  predict_gaussian_latitudes_hemisphere(N,lat.data());
  setup_lat_hemisphere(N,lon,lat.data(),DEG);
}

}
}
}
