/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#include <algorithm>

#include "atlas/meshgen/RGG.hpp"

namespace atlas {
namespace meshgen {
  
T63::T63()
{
  int nlat=32;
  int lon[] = {
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
    125,
    125,
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
double colat[] = {    
    0.03728374374031613592,
    0.08558174883654483256,
    0.13416497894680912140,
    0.18281265245634631400,
    0.23148476959988525947,
    0.28016871368937540376,
    0.32885926587507940955,
    0.37755388050436688996,
    0.42625116887703462520,
    0.47495030929500642314,
    0.52365078451647795177,
    0.57235225266232836461,
    0.62105447864251428491,
    0.66975729540951212115,
    0.71846058092900699776,
    0.76716424390145598888,
    0.81586821458595581991,
    0.86457243871818434400,
    0.91327687336912655169,
    0.96198148405750538714,
    1.01068624269321394316,
    1.05939112608421659445,
    1.10809611483324954584,
    1.15680119250898072458,
    1.20550634501341691340,
    1.25421156009148382360,
    1.30291682694470245529,
    1.35162213592166891019,
    1.40032747826539116787,
    1.44903284590263159437,
    1.49773823126390936977,
    1.54644362712526550752
  };
  
  lat_.resize(2*nlat);
  lon_.resize(2*nlat);
  std::copy( lon, lon+nlat, lon_.begin() );
  std::reverse_copy( lon, lon+nlat, lon_.begin()+nlat );
  std::copy( colat, colat+nlat, lat_.begin() );
  std::reverse_copy( colat, colat+nlat, lat_.begin()+nlat );
  for (int i=0; i<nlat; ++i)
    lat_[i]=M_PI/2.-lat_[i];
  for (int i=nlat; i<2*nlat; ++i)
    lat_[i]=-M_PI/2.+lat_[i];
}

} // namespace meshgen
} // namespace atlas
