/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <sstream>

#include "atlas/atlas_config.h"
#include "atlas/meshgen/RGG.hpp"
#include "atlas/meshgen/EqualAreaPartitioner.hpp"
#include "atlas/Gmsh.hpp"
#include "atlas/BuildEdges.hpp"
#include "atlas/BuildDualMesh.hpp"
#include "atlas/BuildPeriodicBoundaries.hpp"
#include "atlas/testing/UnitTest.hpp"

using namespace atlas;
using namespace atlas::meshgen;

class TestMeshGen: public UnitTest
{
public:
  TestMeshGen(int argc, char *argv[]) : UnitTest(argc,argv) {}
  
  virtual void run_tests()
  {
    test_eq_caps();
    test_partitioner();
    test_rgg_meshgen();
  }
  
  void test_eq_caps()
  {
    std::vector<int>    n_regions;
    std::vector<double> s_cap;

    eq_caps(6, n_regions, s_cap);
    ATLAS_CHECK_EQUAL( n_regions.size(), 3 );
    ATLAS_CHECK_EQUAL( n_regions[0], 1 );
    ATLAS_CHECK_EQUAL( n_regions[1], 4 );
    ATLAS_CHECK_EQUAL( n_regions[2], 1 );

    eq_caps(10, n_regions, s_cap);
    ATLAS_CHECK_EQUAL( n_regions.size(), 4 );
    ATLAS_CHECK_EQUAL( n_regions[0], 1 );
    ATLAS_CHECK_EQUAL( n_regions[1], 4 );
    ATLAS_CHECK_EQUAL( n_regions[2], 4 );
    ATLAS_CHECK_EQUAL( n_regions[3], 1 );
  }
  
  void test_partitioner()
  {    
    // 12 partitions
    {
      EqualAreaPartitioner partitioner(12);
      ATLAS_CHECK_EQUAL( partitioner.nb_bands(),    4 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(0), 1 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(1), 5 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(2), 5 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(3), 1 );
    }
    
    // 24 partitions
    {
      EqualAreaPartitioner partitioner(24);
      ATLAS_CHECK_EQUAL( partitioner.nb_bands(),     5 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(0),  1 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(1),  6 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(2), 10 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(3),  6 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(4),  1 );
    }
    
    // 48 partitions
    {
      EqualAreaPartitioner partitioner(48);
      ATLAS_CHECK_EQUAL( partitioner.nb_bands(),     7 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(0),  1 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(1),  6 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(2), 11 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(3), 12 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(4), 11 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(5),  6 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(6),  1 );
    }
    
    // 96 partitions
    {
      EqualAreaPartitioner partitioner(96);
      ATLAS_CHECK_EQUAL( partitioner.nb_bands(),    10 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(0),  1 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(1),  6 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(2), 11 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(3), 14 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(4), 16 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(5), 16 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(6), 14 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(7), 11 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(8),  6 );
      ATLAS_CHECK_EQUAL( partitioner.nb_regions(9),  1 );
    }
  }
  
  void test_rgg_meshgen()
  {
    RGGMeshGenerator generate;
    generate.options.set("nb_parts",20);

    for( int p=0; p<generate.options.get<int>("nb_parts"); ++p)
    {
      generate.options.set("part",p);
      Mesh* m = generate( T63() );
      std::stringstream filename; filename << "d" << p << ".msh";
      Gmsh::write(*m,filename.str());
    }
  
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
  }
};


int main(int argc, char *argv[])
{
  TestMeshGen tests(argc,argv);
  tests.start();
  return tests.exit_status();
}
