// (C) Copyright 1996-2014 ECMWF.

#include "atlas/atlas_config.h"

#include "atlas/MPL.hpp"
#include "atlas/meshgen/RGG.hpp"
#include "atlas/meshgen/EqualAreaPartitioner.hpp"
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
  
  int N=6;
  std::vector<int> n_regions;
  std::vector<double> s_cap;
  eq_caps(N, n_regions, s_cap);
  for (int n=0; n<n_regions.size(); ++n){
    std::cout << n_regions[n] << std::endl;
  }
  for (int n=0; n<s_cap.size(); ++n){
    std::cout << s_cap[n] << std::endl;
  }
  
  std::vector<double> xmin(N);
  std::vector<double> xmax(N);
  std::vector<double> ymin(N);
  std::vector<double> ymax(N);
  eq_regions(N,xmin.data(),xmax.data(),ymin.data(),ymax.data());
  for (int n=0; n<N; ++n)
  {
    std::cout << n<<" : [ ["<<xmin[n]<<","<<ymin[n]<<"] , ["<<xmax[n]<<","<<ymax[n] <<"] ]" << std::endl;
  }
  
  EqualAreaPartitioner partitioner(1028);
  
  std::cout << "part = " << partitioner.partition(4.5,-0.3) << std::endl;
  
  MPL::finalize();
  return 0;
}
