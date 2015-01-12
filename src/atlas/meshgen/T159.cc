/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/grids/GaussianLatitudes.h"
#include "atlas/meshgen/RGG.h"

namespace atlas {
namespace meshgen {

T159::T159()
{
  int N=80;
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
    200,
    216,
    216,
    216,
    225,
    225,
    240,
    240,
    240,
    243,
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
    288,
    288,
    300,
    300,
    300,
    300,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320,
    320
  };
  std::vector<double> lats(N);
  grids::gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon,DEG);
}

} // namespace meshgen
} // namespace atlas

