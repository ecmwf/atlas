// TL215

#include <eckit/log/Log.h>
#include "atlas/grids/reduced_gg/Grids.h"

namespace atlas {
namespace grids {
namespace reduced_gg {

N108::regist N108_builders;

void N108::construct()
{
  int N=108;
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
   216,
   216,
   216,
   225,
   240,
   240,
   240,
   243,
   250,
   256,
   256,
   270,
   270,
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
   375,
   375,
   375,
   375,
   384,
   384,
   384,
   400,
   400,
   400,
   400,
   400,
   400,
   405,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432,
   432
  };
  std::vector<double> lat(N);
  eckit::Log::warning() << className() << " uses predicted gaussian latitudes" << std::endl;
  predict_gaussian_latitudes_hemisphere(N,lat.data());
  setup_lat_hemisphere(N,lon,lat.data(),DEG);
}

}
}
}
