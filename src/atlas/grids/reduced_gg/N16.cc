// TL31

#include "atlas/grids/reduced_gg/Grids.h"

namespace atlas {
namespace grids {
namespace reduced_gg {

N16::regist N16_builders;

void N16::construct()
{
  int N=16;
  int lon[] = {
      20,
      27,
      32,
      40,
      45,
      48,
      60,
      60,
      64,
      64,
      64,
      64,
      64,
      64,
      64,
      64
  };
  std::vector<double> lat(N);
  eckit::Log::warning() << className() << " uses predicted gaussian latitudes";
  predict_gaussian_latitudes_hemisphere(N,lat.data());
  setup_lat_hemisphere(N,lon,lat.data(),DEG);
}

}
}
}
