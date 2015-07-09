/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grids/rgg/OctahedralRGG.h"

namespace atlas {
namespace grids {
namespace rgg {

//------------------------------------------------------------------------------------------------------

std::vector<long> OctahedralRGG::computePL(const size_t N, const size_t start)
{
  std::vector<long> nlon(N);
  for(size_t jlat=0; jlat < N; ++jlat)
  {
    nlon[jlat] = start + 4*jlat;
  }
  return nlon;
}

OctahedralRGG::OctahedralRGG(const size_t N, const size_t octahedralPoleStart)
{
  construct(N,octahedralPoleStart);
  set_typeinfo();
}

OctahedralRGG::OctahedralRGG( const eckit::Parametrisation& params)
{
    size_t N;
    params.get("N",N);

    size_t octahedralPoleStart = 20;
    if(params.has("OctahedralPoleStart")) {
        params.get("OctahedralPoleStart",octahedralPoleStart);
    }

    construct(N,octahedralPoleStart);
    set_typeinfo();
}

void OctahedralRGG::construct(const size_t N, const size_t start)
{
  std::vector<long> nlon = computePL(N,start);
  setup_N_hemisphere(N,nlon.data());
  ReducedGrid::N_ = nlat()/2;
}

void OctahedralRGG::set_typeinfo()
{
    std::ostringstream s;
    s << "oct.N"<< N();
    shortName_ = s.str();
    grid_type_ = ReducedGaussianGrid::grid_type_str();
}

eckit::ConcreteBuilderT1<Grid,OctahedralRGG> builder_OctahedralRGG (OctahedralRGG::grid_type_str());

//------------------------------------------------------------------------------------------------------

} // namespace rgg
} // namespace grids
} // namespace atlas
