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

  RGGMeshGenerator generate;
  generate.options.set("nb_parts",20);
    
  generate.options.set("part",1);
  Mesh* m1 = generate( T63() );
  Gmsh::write(*m1,"w63_1.msh");

  generate.options.set("part",2);
  Mesh* m2 = generate( T63() );
  Gmsh::write(*m2,"w63_2.msh");

  generate.options.set("part",11);
  Mesh* m11 = generate( T63() );
  Gmsh::write(*m11,"w63_11.msh");
  
  // build_periodic_boundaries(*mesh);
  // build_edges(*mesh);
  // build_dual_mesh(*mesh);
  
  
  // Gmsh::write(*mesh,"w63.msh");
#if 0
  int N=6;
  std::vector<int> n_regions;
  std::vector<double> s_cap;
  eq_caps(N, n_regions, s_cap);
  std::cout << "regions" << std::endl;
  for (int n=0; n<n_regions.size(); ++n){
    std::cout << n_regions[n] << std::endl;
  }
  std::cout << "caps" << std::endl;
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
  
  EqualAreaPartitioner partitioner(6);
  
  double x =  M_PI/4.;
  double y = 0.;
  std::cout << "part = " << partitioner.partition(x,y) << std::endl;
  std::cout << "band = " << partitioner.band(y) << std::endl;
  std::cout << "sector = " << partitioner.sector( partitioner.band(y),x ) <<std::endl;
  
#endif
  MPL::finalize();
  return 0;
}
