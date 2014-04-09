// (C) Copyright 1996-2014 ECMWF.

#include "atlas/atlas_config.h"

#include "atlas/MPL.hpp"
#include "atlas/meshgen/RGG.hpp"
#include "atlas/Gmsh.hpp"
#include "atlas/BuildEdges.hpp"
#include "atlas/BuildDualMesh.hpp"
#include "atlas/BuildPeriodicBoundaries.hpp"

using namespace atlas;
using namespace atlas::meshgen;

int main(int argc, char *argv[])
{
  MPL::init();

  RGGMeshGenerator gen;
  Mesh* mesh = gen.generate( T159() );
  
  // build_periodic_boundaries(*mesh);
  // build_edges(*mesh);
  // build_dual_mesh(*mesh);
  
  
  Gmsh::write(*mesh,"t159.msh");
  
  MPL::finalize();
  return 0;
}
