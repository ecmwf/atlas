// TL31

#include "atlas/grid/predefined/rgg/rgg.h"

namespace atlas {
namespace grid {
namespace predefined {
namespace rgg {

eckit::ConcreteBuilderT1<Grid,N16> builder_N16_deprecated (N16::grid_type_str());
eckit::ConcreteBuilderT1<Grid,N16> builder_N16 ("N16");

void N16::construct()
{
  int N=16;
  long lon[] = {
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
  std::vector<double> lats(N);
  gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon,internals::DEG);
}

} // namespace rgg
} // namespace predefined
} // namespace grid
} // namespace atlas
