// TL95

#include "atlas/grids/rgg/rgg.h"

namespace atlas {
namespace grids {
namespace rgg {

eckit::ConcreteBuilderT1<Grid,N48> builder_N48 (N48::grid_type_str());
eckit::ConcreteBuilderT1<Grid,N48> builder_TL95("rgg.TL95");

void N48::construct()
{
  int N=48;
  int lon[] = {
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
  setup_lat_hemisphere(N,lats.data(),lon,DEG);
}

} // namespace rgg
} // namespace grids
} // namespace atlas
