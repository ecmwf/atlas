// TL47

#include "atlas/grids/rgg/rgg.h"

namespace atlas {
namespace grids {
namespace rgg {

eckit::ConcreteBuilderT1<Grid,N24> builder_N24 (N24::grid_type_str());

void N24::construct()
{
  int N=24;
  long lon[] = {
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
  std::vector<double> lats(N);
  gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon,DEG);
}

} // namespace rgg
} // namespace grids
} // namespace atlas
