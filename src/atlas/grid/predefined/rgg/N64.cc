// TL127

#include "atlas/grid/predefined/rgg/rgg.h"

namespace atlas {
namespace grid {
namespace predefined {
namespace rgg {

eckit::ConcreteBuilderT1<Grid,N64> builder_N64_deprecated (N64::grid_type_str());
eckit::ConcreteBuilderT1<Grid,N64> builder_N64 ("N64");

void N64::construct()
{
  int N=64;
  long lon[] = {
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
  std::vector<double> lats(N);
  gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon,internals::DEG);
}

} // namespace rgg
} // namespace predefined
} // namespace grid
} // namespace atlas
