// TL31

#include "atlas/grids/reduced_gg/Grids.h"

namespace atlas {
namespace grids {
namespace reduced_gg {

void N16::construct()
{
  int N=16;
  int lon[] = {
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
  double lat[] = {
    85.7605871204438159,
    80.2687790722500125,
    74.7445403686357679,
    69.2129761693708616,
    63.6786355610968613,
    58.1429540492032828,
    52.6065260343452650,
    47.0696420596876806,
    41.5324612466560765,
    35.9950784112715994,
    30.4575539611520938,
    24.9199286299486111,
    19.3822313464343878,
    13.8444837343848572,
     8.3067028565188039,
     2.7689030077360099
  };
  setup_lat_hemisphere(N,lon,lat,DEG);
}

} // namespace reduced_gg
} // namespace grids
} // namespace atlas
