// (C) Copyright 1996-2014 ECMWF.

#include "atlas/atlas_config.h"

#include "atlas/Gmsh.hpp"
#include "atlas/BuildEdges.hpp"
#include "atlas/BuildDualMesh.hpp"
#include "atlas/BuildPeriodicBoundaries.hpp"
#include "atlas/Partitioner.hpp"
#include "atlas/MPL.hpp"

using namespace atlas;
int main(int argc, char *argv[])
{
  MPL::init();

  Mesh& mesh = Gmsh::read( "T47.msh" );

  build_periodic_boundaries(mesh);
  build_edges(mesh);
  build_dual_mesh(mesh);

  Gmsh::write(mesh,"bla.msh");
  
  MPL::finalize();
  return 0;
}
