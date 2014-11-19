// TL127

#include "atlas/grids/reduced_gg/Grids.h"

namespace atlas {
namespace grids {
namespace reduced_gg {

N64::regist N64_builders;

void N64::construct()
{
  int N=64;
  int lon[] = {
    20,
    25,
    36,
    40,
    45,
    54,
    60,
    64,
    72,
    75,
    80,
    90,
    96,
   100,
   108,
   120,
   120,
   125,
   135,
   135,
   144,
   150,
   160,
   160,
   180,
   180,
   180,
   180,
   192,
   192,
   200,
   200,
   216,
   216,
   216,
   216,
   225,
   225,
   225,
   240,
   240,
   240,
   240,
   243,
   250,
   250,
   250,
   250,
   256,
   256,
   256,
   256,
   256,
   256,
   256,
   256,
   256,
   256,
   256,
   256,
   256,
   256,
   256,
   256
  };
  std::vector<double> lat(N);
  eckit::Log::warning() << className() << " uses predicted gaussian latitudes" << std::endl;
  predict_gaussian_latitudes_hemisphere(N,lat.data());
  setup_lat_hemisphere(N,lon,lat.data(),DEG);
}

}
}
}
