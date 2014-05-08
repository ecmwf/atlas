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

#include "atlas/MPL.hpp"
#include "atlas/meshgen/RGG.hpp"
#include "atlas/meshgen/EqualAreaPartitioner.hpp"
#include "atlas/Gmsh.hpp"
#include "atlas/BuildEdges.hpp"
#include "atlas/BuildDualMesh.hpp"
#include "atlas/BuildPeriodicBoundaries.hpp"
#include <sstream>
using namespace atlas;
using namespace atlas::meshgen;

int main(int argc, char *argv[])
{
  MPL::init();

  RGGMeshGenerator generate;
  generate.options.set("nb_parts",20);

  for( int p=0; p<generate.options.get<int>("nb_parts"); ++p)
  {
    generate.options.set("part",p);
    Mesh* m = generate( T63() );
    std::stringstream filename;
    filename << "d" << p << ".msh";
    Gmsh::write(*m,filename.str());
  }
  // 
  // generate.options.set("part",0);
  // Mesh* m0 = generate( T63() );
  // Gmsh::write(*m0,"d0.msh");
  //   
  // generate.options.set("part",1);
  // Mesh* m1 = generate( T63() );
  // Gmsh::write(*m1,"d1.msh");
  // 
  // generate.options.set("part",2);
  // Mesh* m2 = generate( T63() );
  // Gmsh::write(*m2,"d2.msh");
  // 
  // generate.options.set("part",3);
  // Mesh* m3 = generate( T63() );
  // Gmsh::write(*m3,"d3.msh");
  // 
  // generate.options.set("part",4);
  // Mesh* m4 = generate( T63() );
  // Gmsh::write(*m4,"d4.msh");
  // 
  // generate.options.set("part",5);
  // Mesh* m5 = generate( T63() );
  // Gmsh::write(*m5,"d5.msh");
  //  
  // generate.options.set("part",6);
  // Mesh* m6 = generate( T63() );
  // Gmsh::write(*m6,"d6.msh");
  //  
  // generate.options.set("part",7);
  // Mesh* m7 = generate( T63() );
  // Gmsh::write(*m7,"d7.msh");
  //  
  // generate.options.set("part",8);
  // Mesh* m8 = generate( T63() );
  // Gmsh::write(*m8,"d8.msh");
  // 
  // generate.options.set("part",9);
  // Mesh* m9 = generate( T63() );
  // Gmsh::write(*m9,"d9.msh");
  // 
  // generate.options.set("part",10);
  // Mesh* m10 = generate( T63() );
  // Gmsh::write(*m10,"d10.msh");
  // 
  // generate.options.set("part",11);
  // Mesh* m11 = generate( T63() );
  // Gmsh::write(*m11,"d11.msh");
  // 
  // generate.options.set("nb_parts",1);
  // generate.options.set("part",0);
  // Mesh* m = generate( T63() );
  // Gmsh::write(*m,"w63.msh");
  
  
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
