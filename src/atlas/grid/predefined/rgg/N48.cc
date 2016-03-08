// TL95

#include "atlas/grid/predefined/rgg/rgg.h"

namespace atlas {
namespace grid {
namespace predefined {
namespace rgg {

eckit::ConcreteBuilderT1<Grid,N48> builder_N48_deprecated (N48::grid_type_str());
eckit::ConcreteBuilderT1<Grid,N48> builder_N48 ("N48");

void N48::construct()
{
  int N=48;
  long lon[] = {
     20,
     25,
     36,
     40,
     45,
     50,
     60,
     60,
     72,
     75,
     80,
     90,
     96,
    100,
    108,
    120,
    120,
    120,
    128,
    135,
    144,
    144,
    160,
    160,
    160,
    160,
    160,
    180,
    180,
    180,
    180,
    180,
    192,
    192,
    192,
    192,
    192,
    192,
    192,
    192,
    192,
    192,
    192,
    192,
    192,
    192,
    192,
    192
  };
  std::vector<double> lats(N);
  gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon);
}

} // namespace rgg
} // namespace predefined
} // namespace grid
} // namespace atlas
