// TL191

#include "atlas/grid/predefined/rgg/rgg.h"

namespace atlas {
namespace grid {
namespace predefined {
namespace rgg {

eckit::ConcreteBuilderT1<Grid,N96> builder_N96 ("N96");
eckit::ConcreteBuilderT1<Grid,N96> builder_N96_deprecated (N96::grid_type_str());

void N96::construct()
{
  int N=96;
  long lon[] = {
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
  std::vector<double> lats(N);
  gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon,internals::DEG);
}

} // namespace rgg
} // namespace predefined
} // namespace grid
} // namespace atlas
