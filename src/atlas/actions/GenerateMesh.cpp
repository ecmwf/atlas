/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mpl/MPL.hpp"
#include "atlas/meshgen/RGG.hpp"
#include "atlas/meshgen/RGGMeshGenerator.hpp"
#include "atlas/actions/GenerateMesh.hpp"
using namespace atlas::meshgen;

namespace atlas {
namespace actions {

namespace {
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

}

Mesh* generate_reduced_gaussian_grid( const std::string& identifier )
{
  RGGMeshGenerator generate;
  generate.options.set( "nb_parts", MPL::size() );
  generate.options.set( "part"    , MPL::rank() );

  Mesh* mesh = 0;

  if     ( identifier == "T63"   ) mesh = generate(T63()  );
  else if( identifier == "T95"   ) mesh = generate(T95()  );
  else if( identifier == "T159"  ) mesh = generate(T159() );
  else if( identifier == "T255"  ) mesh = generate(T255() );
  else if( identifier == "T511"  ) mesh = generate(T511() );
  else if( identifier == "T1279" ) mesh = generate(T1279());
  else if( identifier == "T2047" ) mesh = generate(T2047());
  else if( identifier == "T3999" ) mesh = generate(T3999());
  else if( identifier == "T7999" ) mesh = generate(T7999());
  else if( identifier == "D5")
  {
    int lon[] = {4,6,8,8,8};
    TestGrid grid(5,lon);

    mesh = generate( grid );
  }
  else throw eckit::BadParameter("Cannot generate mesh "+identifier,Here());

  return mesh;
}

// ------------------------------------------------------------------

Mesh* generate_latlon_grid( int nlon, int nlat )
{
  if( nlon%2 != 0 ) throw eckit::BadParameter("nlon must be even number",Here());
  if( nlat%2 != 0 ) throw eckit::BadParameter("nlat must be even number",Here());

  RGGMeshGenerator generate;
  generate.options.set( "nb_parts", MPL::size() );
  generate.options.set( "part"    , MPL::rank() );
  Mesh* mesh = generate(GG(nlon,nlat/2));
  return mesh;
}

// ------------------------------------------------------------------

// C wrapper interfaces to C++ routines
Mesh* atlas__generate_reduced_gaussian_grid (char* identifier)
{
  return generate_reduced_gaussian_grid( std::string(identifier) );
}

Mesh* atlas__generate_latlon_grid ( int nlon, int nlat )
{
  return generate_latlon_grid( nlon, nlat );
}

// ------------------------------------------------------------------


} // namespace actions
} // namespace atlas
