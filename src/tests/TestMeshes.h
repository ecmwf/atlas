/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/atlas_config.h"
#include "atlas/meshgen/RGG.h"
#include "atlas/meshgen/RGGMeshGenerator.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mpl/MPL.h"

using namespace atlas;

namespace atlas {
namespace test {
  
class TestGrid: public meshgen::RGG {
public:
  TestGrid(int nlat, int lon[]);
};

TestGrid::TestGrid(int nlat, int lon[]) : RGG()
{
  /*
  First prediction of colatitudes
  */
  std::vector<double> colat(nlat);
  double z;
  for( int i=0; i<nlat; ++i )
  {
    z = (4.*(i+1.)-1.)*M_PI/(4.*2.*nlat+2.);
    colat[i] = z+1./(tan(z)*(8.*(2.*nlat)*(2.*nlat)));
  }
  /*
  Fill in final structures
  */
  lat_.resize(2*nlat);
  lon_.resize(2*nlat);
  std::copy( lon, lon+nlat, lon_.begin() );
  std::reverse_copy( lon, lon+nlat, lon_.begin()+nlat );
  std::copy( colat.begin(), colat.begin()+nlat, lat_.begin() );
  std::reverse_copy( colat.begin(), colat.begin()+nlat, lat_.begin()+nlat );
  for (int i=0; i<nlat; ++i)
    lat_[i]=M_PI/2.-lat_[i];
  for (int i=nlat; i<2*nlat; ++i)
    lat_[i]=-M_PI/2.+lat_[i];
}

Mesh::Ptr generate_mesh( const meshgen::RGG& rgg )
{
  meshgen::RGGMeshGenerator generator;
  generator.options.set("nb_parts",MPL::size());
  generator.options.set("part",MPL::rank());
  return Mesh::Ptr( generator.generate( rgg ) );
}

Mesh::Ptr generate_mesh(int nlat, int lon[] )
{
  return generate_mesh( TestGrid(nlat,lon) );
}


} // end namespace test
} // end namespace atlas
