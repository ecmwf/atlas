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
#include <algorithm>
#include <iomanip>
#define BOOST_TEST_MODULE TestRGG
#define BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY
#include "ecbuild/boost_test_framework.h"

#include "atlas/mpl/MPL.hpp"
#include "atlas/atlas_config.h"
#include "atlas/meshgen/RGG.hpp"
#include "atlas/meshgen/RGGMeshGenerator.hpp"
#include "atlas/meshgen/EqualAreaPartitioner.hpp"
#include "atlas/io/Gmsh.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Metadata.hpp"
#include "atlas/util/ArrayView.hpp"
#include "atlas/util/IndexView.hpp"
#include "atlas/actions/BuildParallelFields.hpp"
#include "atlas/mesh/Parameters.hpp"
#include "atlas/mesh/Util.hpp"

using namespace atlas;
using namespace atlas::meshgen;

#define DISABLE if(0)
#define ENABLE if(1)

namespace atlas {
namespace test {
  
class DebugMesh:   public RGG { public: DebugMesh(); };
DebugMesh::DebugMesh()
{
  int nlat=5;
  int lon[] = {
    6,
    10,
    18,
    22,
    22,
  };
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

  
class MinimalMesh:   public RGG { public: MinimalMesh(int nlat, int lon[]); };
MinimalMesh::MinimalMesh(int nlat, int lon[])
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

double compute_latlon_area(Mesh& mesh)
{
  FunctionSpace& nodes  = mesh.function_space("nodes");
  FunctionSpace& quads  = mesh.function_space("quads");
  FunctionSpace& triags = mesh.function_space("triags");
  ArrayView<double,2> latlon  ( nodes.field("coordinates") );
  IndexView<int,2> quad_nodes ( quads. field("nodes") );
  IndexView<int,2> triag_nodes( triags.field("nodes") );
  double area=0;
  for( int e=0; e<quads.extents()[0]; ++e )
  {
    int n0 = quad_nodes(e,0);
    int n1 = quad_nodes(e,1);
    int n2 = quad_nodes(e,2);
    int n3 = quad_nodes(e,3);
    double x0=latlon(n0,XX), x1=latlon(n1,XX), x2=latlon(n2,XX), x3=latlon(n3,XX);
    double y0=latlon(n0,YY), y1=latlon(n1,YY), y2=latlon(n2,YY), y3=latlon(n3,YY);
    area += std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
    area += std::abs( x2*(y3-y0)+x3*(y0-y2)+x0*(y2-y3) )*0.5;
  }
  for( int e=0; e<triags.extents()[0]; ++e )
  {
    int n0 = triag_nodes(e,0);
    int n1 = triag_nodes(e,1);
    int n2 = triag_nodes(e,2);
    double x0=latlon(n0,XX), x1=latlon(n1,XX), x2=latlon(n2,XX);
    double y0=latlon(n0,YY), y1=latlon(n1,YY), y2=latlon(n2,YY);
    area += std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
  }
  return area;
}


} // end namespace test
} // end namespace atlas

BOOST_AUTO_TEST_CASE( init ) { MPL::init(); }

BOOST_AUTO_TEST_CASE( test_eq_caps )
{
  std::vector<int>    n_regions;
  std::vector<double> s_cap;

  eq_caps(6, n_regions, s_cap);
  BOOST_CHECK_EQUAL( n_regions.size(), 3 );
  BOOST_CHECK_EQUAL( n_regions[0], 1 );
  BOOST_CHECK_EQUAL( n_regions[1], 4 );
  BOOST_CHECK_EQUAL( n_regions[2], 1 );

  eq_caps(10, n_regions, s_cap);
  BOOST_CHECK_EQUAL( n_regions.size(), 4 );
  BOOST_CHECK_EQUAL( n_regions[0], 1 );
  BOOST_CHECK_EQUAL( n_regions[1], 4 );
  BOOST_CHECK_EQUAL( n_regions[2], 4 );
  BOOST_CHECK_EQUAL( n_regions[3], 1 );
}
  
BOOST_AUTO_TEST_CASE( test_partitioner )
{
  // 12 partitions
  {
    EqualAreaPartitioner partitioner(12);
    BOOST_CHECK_EQUAL( partitioner.nb_bands(),    4 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(0), 1 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(1), 5 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(2), 5 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(3), 1 );
  }

  // 24 partitions
  {
    EqualAreaPartitioner partitioner(24);
    BOOST_CHECK_EQUAL( partitioner.nb_bands(),     5 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(0),  1 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(1),  6 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(2), 10 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(3),  6 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(4),  1 );
  }

  // 48 partitions
  {
    EqualAreaPartitioner partitioner(48);
    BOOST_CHECK_EQUAL( partitioner.nb_bands(),     7 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(0),  1 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(1),  6 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(2), 11 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(3), 12 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(4), 11 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(5),  6 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(6),  1 );
  }

  // 96 partitions
  {
    EqualAreaPartitioner partitioner(96);
    BOOST_CHECK_EQUAL( partitioner.nb_bands(),    10 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(0),  1 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(1),  6 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(2), 11 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(3), 14 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(4), 16 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(5), 16 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(6), 14 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(7), 11 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(8),  6 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(9),  1 );
  }
}
  
BOOST_AUTO_TEST_CASE( test_rgg_meshgen_one_part )
{
  Mesh* m;
  RGGMeshGenerator generate;
  generate.options.set("nb_parts",1);
  generate.options.set("part",    0);

  ENABLE {
    generate.options.set("three_dimensional",true);
    generate.options.set("include_pole",false);
    m = generate( atlas::test::DebugMesh() );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).extents()[0], 156 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).extents()[0], 134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").extents()[0],  32 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<int>("max_glb_idx"), 156 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<int>("nb_owned"),    156 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<int>("max_glb_idx"), 166 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<int>("nb_owned"),    134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<int>("max_glb_idx"), 166 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<int>("nb_owned"),     32 );
    delete m;
  }

  ENABLE {
    generate.options.set("three_dimensional",false);
    generate.options.set("include_pole",false);
    m = generate( atlas::test::DebugMesh() );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).extents()[0], 166 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).extents()[0], 134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").extents()[0],  32 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<int>("max_glb_idx"), 166 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<int>("nb_owned"),    166 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<int>("max_glb_idx"), 166 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<int>("nb_owned"),    134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<int>("max_glb_idx"), 166 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<int>("nb_owned"),     32 );
    delete m;
  }

  ENABLE {
    generate.options.set("three_dimensional",true);
    generate.options.set("include_pole",true);
    m = generate( atlas::test::DebugMesh() );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).extents()[0], 158 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).extents()[0], 134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").extents()[0],  44 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<int>("max_glb_idx"), 158 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<int>("nb_owned"),    158 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<int>("max_glb_idx"), 178 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<int>("nb_owned"),    134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<int>("max_glb_idx"), 178 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<int>("nb_owned"),     44 );
    delete m;
  }

  Mesh* mesh;

  ENABLE {
    generate.options.set("three_dimensional",false);
    generate.options.set("include_pole",false);
    int nlat=2;
    int lon[] = { 4, 6 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->function_space("nodes" ).extents()[0], 24 );
    BOOST_CHECK_EQUAL( mesh->function_space("quads" ).extents()[0], 14 );
    BOOST_CHECK_EQUAL( mesh->function_space("triags").extents()[0],  4 );

    double max_lat = test::MinimalMesh(nlat,lon).lat(0);
    BOOST_CHECK_CLOSE( test::compute_latlon_area(*mesh), 2.*M_PI*2.*max_lat, 1e-8 );
    Gmsh::write(*mesh,"minimal2.msh");
    delete mesh;
  }
  // 3 latitudes
  ENABLE {
    int nlat=3;
    int lon[] = { 4, 6, 8 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->function_space("nodes" ).extents()[0], 42 );
    BOOST_CHECK_EQUAL( mesh->function_space("quads" ).extents()[0], 28 );
    BOOST_CHECK_EQUAL( mesh->function_space("triags").extents()[0],  8 );
    Gmsh::write(*mesh,"minimal3.msh");
    delete mesh;
  }
  // 4 latitudes
  ENABLE {
    int nlat=4;
    int lon[] = { 4, 6, 8, 10 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->function_space("nodes" ).extents()[0], 64 );
    BOOST_CHECK_EQUAL( mesh->function_space("quads" ).extents()[0], 46 );
    BOOST_CHECK_EQUAL( mesh->function_space("triags").extents()[0], 12 );
    Gmsh::write(*mesh,"minimal4.msh");
    delete mesh;
  }
  // 5 latitudes WIP
  ENABLE {
    int nlat=5;
    int lon[] = { 6, 10, 18, 22, 22 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->function_space("nodes" ).extents()[0], 166 );
    BOOST_CHECK_EQUAL( mesh->function_space("quads" ).extents()[0], 134 );
    BOOST_CHECK_EQUAL( mesh->function_space("triags").extents()[0],  32 );
    Gmsh::write(*mesh,"minimal5.msh");
    delete mesh;
  }
}

BOOST_AUTO_TEST_CASE( test_rgg_meshgen_many_parts )
{
  RGGMeshGenerator generate;
  generate.options.set("nb_parts",20);
  generate.options.set("include_pole",false);
  generate.options.set("three_dimensional",false);

          //  Alternative grid for debugging
          //  int nlat=10;
          //  int lon[] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 };
          //  test::MinimalMesh grid(nlat,lon);
  T63 grid;
//  RegularGrid grid(128,64);
  double max_lat = grid.lat(0);
  double check_area = 2.*M_PI*2.*max_lat;
  double area = 0;
  int nodes[]  = {313,332,336,338,334,337,348,359,360,361,360,360,359,370,321,334,338,335,348,315};
  int quads[]  = {243,290,293,294,291,293,310,320,321,321,320,321,320,331,278,291,294,293,305,245};
  int triags[] = { 42, 13, 12, 13, 12, 14,  0,  1,  0,  1,  1,  0,  1,  0, 14, 12, 13, 11, 14, 42};
  int nb_owned = 0;

  std::vector<int> all_owned    ( grid.ngptot()+grid.nlat()+1, -1 );

  for( int p=0; p<generate.options.get<int>("nb_parts"); ++p)
  {
    generate.options.set("part",p);
    Mesh::Ptr m( generate( grid ) );
    BOOST_TEST_CHECKPOINT("generated grid " << p);
//    IndexView<int,1> ridx( m->function_space("nodes").field("remote_idx") );
    ArrayView<int,1> part( m->function_space("nodes").field("partition") );
    ArrayView<int,1> gidx( m->function_space("nodes").field("glb_idx") );

    area += test::compute_latlon_area(*m);
    if( generate.options.get<int>("nb_parts") == 20 )
    {
      BOOST_CHECK_EQUAL( m->function_space("nodes" ).extents()[0], nodes[p]  );
      BOOST_CHECK_EQUAL( m->function_space("quads" ).extents()[0], quads[p]  );
      BOOST_CHECK_EQUAL( m->function_space("triags").extents()[0], triags[p] );
    }
    //    std::stringstream filename; filename << "T63_p" << p << ".msh";
    //    Gmsh::write(*m,filename.str());


    FunctionSpace& nodes = m->function_space("nodes");
    int nb_nodes = nodes.extents()[0];

    // Test if all nodes are connected
    {
      std::vector<int> node_elem_connections( nb_nodes, 0 );

      for( int f=0; f<m->nb_function_spaces(); ++f )
      {
        FunctionSpace& elements = m->function_space(f);
        if( elements.metadata<int>("type") == Entity::ELEMS )
        {
          int nb_elems = elements.extents()[0];
          IndexView<int,2> elem_nodes ( elements.field("nodes") );
          int nb_nodes_per_elem = elem_nodes.extents()[1];
          for( int jelem=0; jelem<nb_elems; ++jelem )
          {
            for( int jnode=0; jnode<nb_nodes_per_elem; ++jnode )
              node_elem_connections[ elem_nodes(jelem,jnode) ]++;
          }
        }
      }
      for( int jnode=0; jnode<nb_nodes; ++jnode )
      {
        BOOST_CHECK_GT( node_elem_connections[jnode] , 0 );
      }
    }

    // Test if all nodes are owned
    ArrayView<int,1> glb_idx( nodes.field("glb_idx") );
    for( int n=0; n<nb_nodes; ++n )
    {
      int owned = (part(n)==p);
      if( owned )
      {
        ++nb_owned;
        BOOST_CHECK( all_owned[ gidx(n) ] == -1 );
        if( all_owned[ gidx(n)] != -1 )
          std::cout << "node " << gidx(n) << " already visited for" << std::endl;
        all_owned[ gidx(n) ] = part(n);
      }
    }
  }

  for( int gid=1; gid<all_owned.size(); ++gid )
  {
    if( all_owned[gid] == -1 )
    {
      BOOST_ERROR( "node " << gid << " is not owned by anyone" );
    }
//    BOOST_CHECK( all_owned[gid] != -1 );
  }
  BOOST_CHECK_EQUAL( nb_owned, grid.ngptot()+grid.nlat() );

  BOOST_CHECK_CLOSE( area, check_area, 1e-10 );

}

BOOST_AUTO_TEST_CASE( finalize ) { MPL::finalize(); }
