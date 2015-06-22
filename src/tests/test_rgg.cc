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
#include "ecbuild/boost_test_framework.h"

#include "eckit/config/ResourceMgr.h"
#include "eckit/geometry/Point3.h"
#include "atlas/mpi/mpi.h"
#include "atlas/atlas_config.h"
#include "atlas/grids/GaussianLatitudes.h"
#include "atlas/grids/grids.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/meshgen/EqualRegionsPartitioner.h"
#include "atlas/io/Gmsh.h"
#include "atlas/Mesh.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/Metadata.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"
#include "atlas/actions/BuildParallelFields.h"
#include "atlas/Parameters.h"
#include "atlas/grids/rgg/rgg.h"

using namespace atlas;
using namespace atlas::io;
using namespace atlas::meshgen;

#define DISABLE if(0)
#define ENABLE if(1)

namespace atlas {
namespace test {


typedef Metadata Parameters;

class DebugMesh:   public grids::ReducedGaussianGrid { public: DebugMesh(); };
DebugMesh::DebugMesh()
{
  int N=5;
  long lon[] = {
    6,
    10,
    18,
    22,
    22,
  };
  std::vector<double> lat(N);
  grids::gaussian_latitudes_npole_equator(N,lat.data());
  setup_lat_hemisphere(N,lat.data(),lon,DEG);
}


class MinimalMesh:   public grids::ReducedGaussianGrid { public: MinimalMesh(int N, long lon[]); };
MinimalMesh::MinimalMesh(int N, long lon[])
{
  std::vector<double> lat(N);
  grids::gaussian_latitudes_npole_equator(N,lat.data());
  setup_lat_hemisphere(N,lat.data(),lon,DEG);
}

double compute_lonlat_area(Mesh& mesh)
{
  FunctionSpace& nodes  = mesh.function_space("nodes");
  FunctionSpace& quads  = mesh.function_space("quads");
  FunctionSpace& triags = mesh.function_space("triags");
  ArrayView<double,2> lonlat  ( nodes.field("lonlat") );
  IndexView<int,2> quad_nodes ( quads. field("nodes") );
  IndexView<int,2> triag_nodes( triags.field("nodes") );
  double area=0;
  for( int e=0; e<quads.shape(0); ++e )
  {
    int n0 = quad_nodes(e,0);
    int n1 = quad_nodes(e,1);
    int n2 = quad_nodes(e,2);
    int n3 = quad_nodes(e,3);
    double x0=lonlat(n0,LON), x1=lonlat(n1,LON), x2=lonlat(n2,LON), x3=lonlat(n3,LON);
    double y0=lonlat(n0,LAT), y1=lonlat(n1,LAT), y2=lonlat(n2,LAT), y3=lonlat(n3,LAT);
    area += std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
    area += std::abs( x2*(y3-y0)+x3*(y0-y2)+x0*(y2-y3) )*0.5;
  }
  for( int e=0; e<triags.shape(0); ++e )
  {
    int n0 = triag_nodes(e,0);
    int n1 = triag_nodes(e,1);
    int n2 = triag_nodes(e,2);
    double x0=lonlat(n0,LON), x1=lonlat(n1,LON), x2=lonlat(n2,LON);
    double y0=lonlat(n0,LAT), y1=lonlat(n1,LAT), y2=lonlat(n2,LAT);
    area += std::abs( x0*(y1-y2)+x1*(y2-y0)+x2*(y0-y1) )*0.5;
  }
  return area;
}


BOOST_AUTO_TEST_CASE( init ) { eckit::mpi::init(); }

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
  grids::ReducedGrid g;

  // 12 partitions
  {
    EqualRegionsPartitioner partitioner(g,12);
    BOOST_CHECK_EQUAL( partitioner.nb_bands(),    4 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(0), 1 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(1), 5 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(2), 5 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(3), 1 );
  }

  // 24 partitions
  {
    EqualRegionsPartitioner partitioner(g,24);
    BOOST_CHECK_EQUAL( partitioner.nb_bands(),     5 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(0),  1 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(1),  6 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(2), 10 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(3),  6 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(4),  1 );
  }

  // 48 partitions
  {
    EqualRegionsPartitioner partitioner(g,48);
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
    EqualRegionsPartitioner partitioner(g,96);
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
  MeshGenerator::Parameters opts;
  opts.set("nb_parts",1);
  opts.set("part",0);
  ReducedGridMeshGenerator generate(opts);
//  generate.options.set("nb_parts",1);
//  generate.options.set("part",    0);
DISABLE{  // This is all valid for meshes generated with MINIMAL NB TRIAGS
  ENABLE {
    generate.options.set("three_dimensional",true);
    generate.options.set("include_pole",false);
    m = generate( atlas::test::DebugMesh() );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).shape(0), 156 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).shape(0), 134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").shape(0),  32 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<size_t>("max_glb_idx"), 156 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<size_t>("nb_owned"),    156 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("max_glb_idx"), 166 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("nb_owned"),    134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("max_glb_idx"), 166 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("nb_owned"),     32 );
    delete m;
  }

  ENABLE {
    generate.options.set("three_dimensional",false);
    generate.options.set("include_pole",false);
    m = generate( atlas::test::DebugMesh() );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).shape(0), 166 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).shape(0), 134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").shape(0),  32 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<size_t>("max_glb_idx"), 166 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<size_t>("nb_owned"),    166 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("max_glb_idx"), 166 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("nb_owned"),    134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("max_glb_idx"), 166 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("nb_owned"),     32 );
    delete m;
  }

  ENABLE {
    generate.options.set("three_dimensional",true);
    generate.options.set("include_pole",true);
    m = generate( atlas::test::DebugMesh() );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).shape(0), 158 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).shape(0), 134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").shape(0),  44 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<size_t>("max_glb_idx"), 158 );
    BOOST_CHECK_EQUAL( m->function_space("nodes" ).metadata().get<size_t>("nb_owned"),    158 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("max_glb_idx"), 178 );
    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("nb_owned"),    134 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("max_glb_idx"), 178 );
    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("nb_owned"),     44 );
    delete m;
  }

  Mesh* mesh;

  ENABLE {
    generate.options.set("three_dimensional",false);
    generate.options.set("include_pole",false);
    int nlat=2;
    long lon[] = { 4, 6 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->function_space("nodes" ).shape(0), 24 );
    BOOST_CHECK_EQUAL( mesh->function_space("quads" ).shape(0), 14 );
    BOOST_CHECK_EQUAL( mesh->function_space("triags").shape(0),  4 );

    double max_lat = test::MinimalMesh(nlat,lon).lat(0);
    BOOST_CHECK_CLOSE( test::compute_lonlat_area(*mesh), 2.*M_PI*2.*max_lat, 1e-8 );
    Gmsh().write(*mesh,"minimal2.msh");
    delete mesh;
  }
  // 3 latitudes
  ENABLE {
    int nlat=3;
    long lon[] = { 4, 6, 8 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->function_space("nodes" ).shape(0), 42 );
    BOOST_CHECK_EQUAL( mesh->function_space("quads" ).shape(0), 28 );
    BOOST_CHECK_EQUAL( mesh->function_space("triags").shape(0),  8 );
    Gmsh().write(*mesh,"minimal3.msh");
    delete mesh;
  }
  // 4 latitudes
  ENABLE {
    int nlat=4;
    long lon[] = { 4, 6, 8, 10 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->function_space("nodes" ).shape(0), 64 );
    BOOST_CHECK_EQUAL( mesh->function_space("quads" ).shape(0), 46 );
    BOOST_CHECK_EQUAL( mesh->function_space("triags").shape(0), 12 );
    Gmsh().write(*mesh,"minimal4.msh");
    delete mesh;
  }
  // 5 latitudes WIP
  ENABLE {
    int nlat=5;
    long lon[] = { 6, 10, 18, 22, 22 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->function_space("nodes" ).shape(0), 166 );
    BOOST_CHECK_EQUAL( mesh->function_space("quads" ).shape(0), 134 );
    BOOST_CHECK_EQUAL( mesh->function_space("triags").shape(0),  32 );
    Gmsh().write(*mesh,"minimal5.msh");
    delete mesh;
  }
}
}

BOOST_AUTO_TEST_CASE( test_rgg_meshgen_many_parts )
{
  eckit::ResourceMgr::instance().set("atlas.meshgen.partitioner","eqregion");
  ReducedGridMeshGenerator generate;
  generate.options.set("nb_parts",20);
  generate.options.set("include_pole",false);
  generate.options.set("three_dimensional",false);

          //  Alternative grid for debugging
          //  int nlat=10;
          //  long lon[] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 };
          //  test::MinimalMesh grid(nlat,lon);
  grids::rgg::N32 grid;
//  RegularGrid grid(128,64);
  double max_lat = grid.lat(0);
  double check_area = 360.*2.*max_lat;
  double area = 0;
  int nodes[]  = {313,332,336,338,334,337,348,359,360,361,360,360,359,370,321,334,338,335,348,315};
  int quads[]  = {243,290,293,294,291,293,310,320,321,321,320,321,320,331,278,291,294,293,305,245};
  int triags[] = { 42, 13, 12, 13, 12, 14,  0,  1,  0,  1,  1,  0,  1,  0, 14, 12, 13, 11, 14, 42};
  int nb_owned = 0;

  std::vector<int> all_owned    ( grid.npts()+grid.nlat()+1, -1 );

  for( int p=0; p<generate.options.get<size_t>("nb_parts"); ++p)
  {
    DEBUG_VAR(p);
    generate.options.set("part",p);
    Mesh::Ptr m( generate( grid ) );
    DEBUG();
    m->metadata().set("part",p);
    BOOST_TEST_CHECKPOINT("generated grid " << p);
    ArrayView<int,1> part( m->function_space("nodes").field("partition") );
    ArrayView<gidx_t,1> gidx( m->function_space("nodes").field("glb_idx") );

    area += test::compute_lonlat_area(*m);
    DEBUG();

    DISABLE {  // This is all valid for meshes generated with MINIMAL NB TRIAGS
    if( generate.options.get<size_t>("nb_parts") == 20 )
    {
      BOOST_CHECK_EQUAL( m->function_space("nodes" ).shape(0), nodes[p]  );
      BOOST_CHECK_EQUAL( m->function_space("quads" ).shape(0), quads[p]  );
      BOOST_CHECK_EQUAL( m->function_space("triags").shape(0), triags[p] );
    }
    }
    DEBUG();

        std::stringstream filename; filename << "T63.msh";
        Gmsh().write(*m,filename.str());


    FunctionSpace& nodes = m->function_space("nodes");
    int nb_nodes = nodes.shape(0);

    // Test if all nodes are connected
    {
      std::vector<int> node_elem_connections( nb_nodes, 0 );

      for( int f=0; f<m->nb_function_spaces(); ++f )
      {
        FunctionSpace& elements = m->function_space(f);
        if( elements.metadata().get<long>("type") == Entity::ELEMS )
        {
          int nb_elems = elements.shape(0);
          IndexView<int,2> elem_nodes ( elements.field("nodes") );
          int nb_nodes_per_elem = elem_nodes.shape(1);
          for( int jelem=0; jelem<nb_elems; ++jelem )
          {
            for( int jnode=0; jnode<nb_nodes_per_elem; ++jnode )
              node_elem_connections[ elem_nodes(jelem,jnode) ]++;
          }
        }
      }
      for( int jnode=0; jnode<nb_nodes; ++jnode )
      {
        if( node_elem_connections[jnode] == 0 )
        {
          std::stringstream ss; ss << "part " << p << ": node_gid " << gidx(jnode) << " is not connected to any element.";
DISABLE{
          BOOST_ERROR( ss.str() );
}
        }
      }
    }

    // Test if all nodes are owned
    ArrayView<gidx_t,1> glb_idx( nodes.field("glb_idx") );
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
  }
  BOOST_CHECK_EQUAL( nb_owned, grid.npts()+grid.nlat() );

  BOOST_CHECK_CLOSE( area, check_area, 1e-10 );

}


BOOST_AUTO_TEST_CASE( test_reduced_lonlat )
{
  int N=11;
  long lon[] = {
    2,  //90
    6,  //72
    12, //54
    18, //36
    24, //18
    24, //0
    26, //-18
    24, //-36
    0,
    0,
    0
  };
  double lat[] ={
     90,
     72,
     54,
     36,
     18,
     0,
    -18,
    -36,
    -54,
    -72,
    -90
  };
  grids::ReducedGrid grid(N,lat,lon);
  ReducedGridMeshGenerator generate;

  bool three_dimensional = true;

  generate.options.set("three_dimensional",three_dimensional);
  generate.options.set("triangulate",false);

  Mesh::Ptr m (generate(grid));

  ArrayView<double,2> lonlat( m->function_space("nodes").field("lonlat") );
  ArrayView<double,2> xyz( m->function_space("nodes").create_field<double>("xyz",3,IF_EXISTS_RETURN) );
  for( int jnode=0; jnode<lonlat.shape(0); ++jnode )
  {
    eckit::geometry::lonlat_to_3d( lonlat[jnode].data(), xyz[jnode].data() );
  }

  io::Gmsh gmsh;
  if(three_dimensional)
    gmsh.options.set("nodes",std::string("xyz"));
  gmsh.write(*m,"rll.msh");

}

BOOST_AUTO_TEST_CASE( finalize ) { eckit::mpi::finalize(); }

} // namespace test
} // namespace atlas
