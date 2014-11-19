// TL191

#include "atlas/grids/reduced_gg/Grids.h"

namespace atlas {
namespace grids {
namespace reduced_gg {

N96::regist N96_builders;

void N96::construct()
{
  int N=96;
  int lon[] = {
    18,
    25,
    36,
    40,
    45,
    50,
    60,
    64,
    72,
    72,
    80,
    90,
    96,
   100,
   108,
   120,
   120,
   125,
   135,
   144,
   144,
   150,
   160,
   160,
   180,
   180,
   180,
   192,
   192,
   200,
   200,
   216,
   216,
   225,
   225,
   240,
   240,
   240,
   250,
   250,
   256,
   270,
   270,
   270,
   288,
   288,
   288,
   288,
   300,
   300,
   300,
   320,
   320,
   320,
   320,
   320,
   324,
   360,
   360,
   360,
   360,
   360,
   360,
   360,
   360,
   360,
   360,
   360,
   375,
   375,
   375,
   375,
   375,
   375,
   375,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384,
   384
  };
  std::vector<double> lat(N);
  eckit::Log::warning() << className() << " uses predicted gaussian latitudes" << std::endl;
  predict_gaussian_latitudes_hemisphere(N,lat.data());
  setup_lat_hemisphere(N,lon,lat.data(),DEG);
}

}
}
}
