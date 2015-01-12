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

T95::T95()
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
  grids::gaussian_latitudes_npole_equator(N,lats.data());
  setup_lat_hemisphere(N,lats.data(),lon,DEG);
}

} // namespace meshgen
} // namespace atlas
