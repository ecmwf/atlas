/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "atlas/projection.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {

//-----------------------------------------------------------------------------


atlas::Projection::Jacobian
getJacobian (const atlas::Projection & proj, const atlas::PointLonLat & lonlat0)
{
  double dlon = 0.0001, dlat = 0.0001;

  double xy0[2] = {lonlat0.lon () + 0.00, lonlat0.lat () + 0.00};
  double xy1[2] = {lonlat0.lon () + dlon, lonlat0.lat () + 0.00};
  double xy2[2] = {lonlat0.lon () + 0.00, lonlat0.lat () + dlat};

  proj.lonlat2xy (xy0);
  proj.lonlat2xy (xy1);
  proj.lonlat2xy (xy2);

  atlas::Projection::Jacobian jac;

  jac[0] = {(xy1[0] - xy0[0]) / dlon, (xy2[0] - xy0[0]) / dlat};
  jac[1] = {(xy1[1] - xy0[1]) / dlon, (xy2[1] - xy0[1]) / dlat};

  return jac;
}

void printJacobian (const atlas::Projection::Jacobian & jac, const std::string & name)
{
  printf (" %s = \n", name.c_str ());
  printf ("  %12.4f, %12.4f | %12.4f\n", jac[0][0], jac[0][1], sqrt (jac[0][0] * jac[0][0] + jac[0][1] * jac[0][1]));
  printf ("  %12.4f, %12.4f | %12.4f\n", jac[1][0], jac[1][1], sqrt (jac[1][0] * jac[1][0] + jac[1][1] * jac[1][1]));
  printf ("  %12.4f, %12.4f\n", sqrt (jac[0][0] * jac[0][0] + jac[1][0] * jac[1][0]),
                                sqrt (jac[0][1] * jac[0][1] + jac[1][1] * jac[1][1]));
}

template <typename SKIP>
void doTaylorTest (const atlas::StructuredGrid & grid, double dmax, SKIP skip)
{
  const auto & proj = grid.projection ();

  for (int j = 0, jglo = 0; j < grid.ny ( ); j++)
  for (int i = 0; i < grid.nx (j); i++, jglo++)
  if (! skip (grid, i, j))
    {
      double lonlat[2];
      grid.lonlat (i, j, lonlat);

      auto jacA = proj.getJacobianAtLonLat (atlas::PointLonLat (lonlat[0], lonlat[1]));
      auto jacB = getJacobian (proj, atlas::PointLonLat (lonlat[0], lonlat[1]));
      auto jacC = jacB * jacA.inverse () - atlas::Projection::Jacobian::Id ();

      double diff = jacC.norm ();

      if (diff > dmax)
        {
          printf ("------------------\n");
          printf (" i, j = %8d, %8d\n", i, j);
          printJacobian (jacA, "jacA");
          printf ("\n");
          printJacobian (jacB, "jacB");
          printf ("\n");
          printf ("%12.4e\n", diff);
          printf ("\n");
        }

      EXPECT (diff < dmax);
    }

}


void doTaylorTest (const atlas::StructuredGrid & grid, double dmax)
{
  auto noskip = [] (const atlas::StructuredGrid & grid, int i, int j) { return false; };
  doTaylorTest (grid, dmax, noskip);
}

CASE ("test_rotated_schmidt")
{

  const std::vector<int> nx =
  { 
    20, 27, 32, 40, 45, 48, 60, 60, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 60, 60, 48, 45, 40, 32, 27, 20
  };
  
  const double stretch = 2.4, centre[2] = {2.0, 46.7};
  
  std::vector<atlas::grid::Spacing> spacings (nx.size ());
  
  for (int i = 0; i < nx.size (); i++)
    {   
      double lonmax = 360.0 * double (nx[i] - 1) / double (nx[i]);
      spacings[i] = atlas::grid::Spacing (atlas::util::Config ("type", "linear") | atlas::util::Config ("N", nx[i])
                                        | atlas::util::Config ("start", 0) | atlas::util::Config ("end", lonmax));
    }
  
  atlas::StructuredGrid::XSpace xspace (spacings);
  atlas::StructuredGrid::YSpace yspace (atlas::util::Config ("type", "gaussian") | atlas::util::Config ("N", nx.size ()));
  
  auto
  proj = atlas::Projection (atlas::util::Config ("type", "rotated_schmidt") 
                          | atlas::util::Config ("stretching_factor", stretch)
                          | atlas::util::Config ("rotation_angle", 0.0)
                          | atlas::util::Config ("north_pole", std::vector<double>{centre[0], centre[1]})
                          | atlas::util::Config ("arpege", true));
  

  doTaylorTest (atlas::StructuredGrid (xspace, yspace, proj, atlas::Domain (atlas::util::Config ("type", "global"))), 1e-3,
  // Avoid large jumps of pseudo-longitude (ie 359.5 -> 0)
                [] (const atlas::StructuredGrid & grid, int i, int j) { return (i == 0) || (i == grid.nx (j) - 1) || (i == grid.nx (j) / 2); });
}


CASE ("test_lambert_conformal_conic")
{
  const int Nx = 64, Ny = 64, Nux = 53, Nuy = 53; 
  const double DxInMetres = 50000., DyInMetres = 50000.;
  const double LaDInDegrees = 46.2, Latin1InDegrees = 46.2, Latin2InDegrees = 46.2, LoVInDegrees = 2.0;
  const double XMinInMeters = -Nux / 2 * DxInMetres, YMinInMetres = -Nuy / 2 * DyInMetres;
  
  std::vector<atlas::grid::Spacing> spacings (Ny);
  
  for (int i = 0; i < Ny; i++)
    spacings[i] = atlas::grid::Spacing (atlas::util::Config ("type", "linear") | atlas::util::Config ("N", Nx) 
                                      | atlas::util::Config ("start", XMinInMeters) 
                                      | atlas::util::Config ("end", XMinInMeters + (Nx - 1) * DxInMetres));
  
  atlas::StructuredGrid::XSpace xspace (spacings);
  atlas::StructuredGrid::YSpace yspace (atlas::util::Config ("type", "linear") | atlas::util::Config ("N", Ny) 
                                      | atlas::util::Config ("start", YMinInMetres) 
                                      | atlas::util::Config ("end", YMinInMetres + (Ny - 1 ) * DyInMetres));

  auto proj = atlas::Projection (atlas::util::Config ("type", "lambert_conformal_conic") | atlas::util::Config ("longitude0", LoVInDegrees)
                               | atlas::util::Config ("latitude0", LaDInDegrees) | atlas::util::Config ("latitude1", Latin1InDegrees) 
                               | atlas::util::Config ("latitude2", Latin2InDegrees));
  
  doTaylorTest (atlas::StructuredGrid (xspace, yspace, proj, atlas::Domain ()), 1e-6);
}

CASE ("test_lonlat")
{
  doTaylorTest (atlas::StructuredGrid ("N16"), 1e-9);
}



//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
