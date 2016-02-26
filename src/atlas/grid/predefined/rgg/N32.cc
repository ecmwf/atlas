// TL63

#include "atlas/grid/predefined/rgg/rgg.h"

namespace atlas {
namespace grid {
namespace predefined {
namespace rgg {

eckit::ConcreteBuilderT1<Grid,N32> builder_N32_deprecated (N32::grid_type_str());
eckit::ConcreteBuilderT1<Grid,N32> builder_N32 ("N32");

void N32::construct()
{
  int N=32;
  long lon[] = {
     20,
     27,
     36,
     40,
     45,
     50,
     60,
     64,
     72,
     75,
     80,
     90,
     90,
     96,
    100,
    108,
    108,
    120,
    120,
    120,
    128,
    128,
    128,
    128,
    128,
    128,
    128,
    128,
    128,
    128,
    128,
    128
  };
  std::vector<double> lats(N);
  gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon,internals::DEG);
}

} // namespace rgg
} // namespace predefined
} // namespace grid
} // namespace atlas
