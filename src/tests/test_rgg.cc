/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/atlas.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/internals/atlas_config.h"
#include "atlas/grid/GaussianLatitudes.h"
#include "atlas/grid/grids.h"
#include "atlas/mesh/generators/ReducedGridMeshGenerator.h"
#include "atlas/grid/partitioners/EqualRegionsPartitioner.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/util/Config.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Elements.h"
#include "atlas/field/Field.h"
#include "atlas/util/Metadata.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/internals/Parameters.h"
#include "atlas/util/Config.h"
#include "atlas/grid/predefined/rgg/rgg.h"
#include "atlas/internals/Bitflags.h"
#include "atlas/internals/Debug.h"

using namespace eckit;
using namespace atlas;
using namespace atlas::internals;
using namespace atlas::grid::partitioners;
using namespace atlas::mesh::generators;
using namespace atlas::util::io;


#define DISABLE if(0)
#define ENABLE if(1)

namespace atlas {
namespace test {

class DebugMesh:   public grid::ReducedGaussianGrid { public: DebugMesh(); };
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
  grid::gaussian_latitudes_npole_equator(N,lat.data());
  setup_lat_hemisphere(N,lat.data(),lon,internals::DEG);
}


class MinimalMesh:   public grid::ReducedGaussianGrid { public: MinimalMesh(int N, long lon[]); };
MinimalMesh::MinimalMesh(int N, long lon[])
{
  std::vector<double> lat(N);
  grid::gaussian_latitudes_npole_equator(N,lat.data());
  setup_lat_hemisphere(N,lat.data(),lon,internals::DEG);
}

double compute_lonlat_area(mesh::Mesh& mesh)
{
  mesh::Nodes& nodes  = mesh.nodes();
  mesh::Elements& quads  = mesh.cells().elements(0);
  mesh::Elements& triags = mesh.cells().elements(1);
  array::ArrayView<double,2> lonlat  ( nodes.lonlat() );

  const mesh::Elements::Connectivity& quad_nodes  = quads.node_connectivity();
  const mesh::Elements::Connectivity& triag_nodes = triags.node_connectivity();

  double area=0;
  for(size_t e = 0; e < quads.size(); ++e)
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
  for(size_t e = 0; e < triags.size(); ++e)
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

struct GlobalFixture {
    GlobalFixture()  { atlas_init(boost::unit_test::framework::master_test_suite().argc,
                                  boost::unit_test::framework::master_test_suite().argv); }
    ~GlobalFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( GlobalFixture );

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
  grid::ReducedGrid g;

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
  mesh::Mesh* m;
  util::Config opts;
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
    BOOST_CHECK_EQUAL( m->nodes().size(), 156 );
    BOOST_CHECK_EQUAL( m->cells().elements(0).size(), 134 );
    BOOST_CHECK_EQUAL( m->cells().elements(1).size(),  32 );
//    BOOST_CHECK_EQUAL( m->nodes().metadata().get<size_t>("nb_owned"),    156 );
//    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("nb_owned"),    134 );
//    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("nb_owned"),     32 );
    delete m;
  }

  ENABLE {
    generate.options.set("three_dimensional",false);
    generate.options.set("include_pole",false);
    m = generate( atlas::test::DebugMesh() );
    BOOST_CHECK_EQUAL( m->nodes().size(), 166 );
    BOOST_CHECK_EQUAL( m->cells().elements(0).size(), 134 );
    BOOST_CHECK_EQUAL( m->cells().elements(1).size(),  32 );
//    BOOST_CHECK_EQUAL( m->nodes().metadata().get<size_t>("nb_owned"),    166 );
//    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("nb_owned"),    134 );
//    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("nb_owned"),     32 );
    delete m;
  }

  ENABLE {
    generate.options.set("three_dimensional",true);
    generate.options.set("include_pole",true);
    m = generate( atlas::test::DebugMesh() );
    BOOST_CHECK_EQUAL( m->nodes().size(), 158 );
    BOOST_CHECK_EQUAL( m->cells().elements(0).size(), 134 );
    BOOST_CHECK_EQUAL( m->cells().elements(1).size(),  44 );
//    BOOST_CHECK_EQUAL( m->nodes().metadata().get<size_t>("nb_owned"),    158 );
//    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("nb_owned"),    134 );
//    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("nb_owned"),     44 );
    delete m;
  }

  mesh::Mesh* mesh;

  ENABLE {
    generate.options.set("three_dimensional",false);
    generate.options.set("include_pole",false);
    int nlat=2;
    long lon[] = { 4, 6 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->nodes().size(), 24 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(0).size(), 14 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(1).size(),  4 );

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
    BOOST_CHECK_EQUAL( mesh->nodes().size(), 42 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(0).size(), 28 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(1).size(),  8 );
    Gmsh().write(*mesh,"minimal3.msh");
    delete mesh;
  }
  // 4 latitudes
  ENABLE {
    int nlat=4;
    long lon[] = { 4, 6, 8, 10 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->nodes().size(), 64 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(0).size(), 46 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(1).size(), 12 );
    Gmsh().write(*mesh,"minimal4.msh");
    delete mesh;
  }
  // 5 latitudes WIP
  ENABLE {
    int nlat=5;
    long lon[] = { 6, 10, 18, 22, 22 };
    mesh = generate( test::MinimalMesh(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->nodes().size(), 166 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(0).size(), 134 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(1).size(),  32 );
    Gmsh().write(*mesh,"minimal5.msh");
    delete mesh;
  }
}
}

BOOST_AUTO_TEST_CASE( test_rgg_meshgen_many_parts )
{
  BOOST_CHECK( PartitionerFactory::has("EqualRegions") );
  eckit::ResourceMgr::instance().set("atlas.meshgen.partitioner","EqualRegions");
  ReducedGridMeshGenerator generate;
  generate.options.set("nb_parts",20);
  generate.options.set("include_pole",false);
  generate.options.set("three_dimensional",false);

          //  Alternative grid for debugging
          //  int nlat=10;
          //  long lon[] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 };
          //  test::MinimalMesh grid(nlat,lon);
  grid::predefined::rgg::N32 grid;
//  RegularGrid grid(128,64);
  double max_lat = grid.lat(0);
  double check_area = 360.*2.*max_lat;
  double area = 0;
  int nodes[]  = {313,332,336,338,334,337,348,359,360,361,360,360,359,370,321,334,338,335,348,315};
  int quads[]  = {243,290,293,294,291,293,310,320,321,321,320,321,320,331,278,291,294,293,305,245};
  int triags[] = { 42, 13, 12, 13, 12, 14,  0,  1,  0,  1,  1,  0,  1,  0, 14, 12, 13, 11, 14, 42};
  int nb_owned = 0;

  std::vector<int> all_owned    ( grid.npts()+grid.nlat()+1, -1 );

  for(size_t p = 0; p < generate.options.get<size_t>("nb_parts"); ++p)
  {
    DEBUG_VAR(p);
    generate.options.set("part",p);
    mesh::Mesh::Ptr m( generate( grid ) );
    DEBUG();
    m->metadata().set("part",p);
    BOOST_TEST_CHECKPOINT("generated grid " << p);
    array::ArrayView<int,1> part( m->nodes().partition() );
    array::ArrayView<gidx_t,1> gidx( m->nodes().global_index() );

    area += test::compute_lonlat_area(*m);
    DEBUG();

    DISABLE {  // This is all valid for meshes generated with MINIMAL NB TRIAGS
    if( generate.options.get<size_t>("nb_parts") == 20 )
    {
      BOOST_CHECK_EQUAL( m->nodes().size(), nodes[p]  );
      BOOST_CHECK_EQUAL( m->cells().elements(0).size(), quads[p]  );
      BOOST_CHECK_EQUAL( m->cells().elements(1).size(), triags[p] );
    }
    }
    DEBUG();

        std::stringstream filename; filename << "T63.msh";
        Gmsh().write(*m,filename.str());


    mesh::Nodes& nodes = m->nodes();
    int nb_nodes = nodes.size();

    // Test if all nodes are connected
    {
      std::vector<int> node_elem_connections( nb_nodes, 0 );

      const mesh::HybridElements::Connectivity& cell_node_connectivity = m->cells().node_connectivity();
      for( size_t jelem=0; jelem<m->cells().size(); ++jelem )
      {
        for( int jnode=0; jnode<cell_node_connectivity.cols(jelem); ++jnode )
          node_elem_connections[ cell_node_connectivity(jelem,jnode) ]++;
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
    array::ArrayView<gidx_t,1> glb_idx( nodes.global_index() );
    for( int n=0; n<nb_nodes; ++n )
    {
      if( size_t(part(n)) == p )
      {
        ++nb_owned;
        BOOST_CHECK( all_owned[ gidx(n) ] == -1 );
        if( all_owned[ gidx(n)] != -1 )
          std::cout << "node " << gidx(n) << " already visited for" << std::endl;
        all_owned[ gidx(n) ] = part(n);
      }
    }
  }

  for(size_t gid = 1; gid < all_owned.size(); ++gid)
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
  grid::ReducedGrid grid(N,lat,lon);
  ReducedGridMeshGenerator generate;

  bool three_dimensional = true;

  generate.options.set("three_dimensional",three_dimensional);
  generate.options.set("triangulate",false);

  mesh::Mesh::Ptr m (generate(grid));

  mesh::actions::BuildXYZField build_xyz_field("xyz");
  build_xyz_field(*m);

  util::io::Gmsh gmsh;
  if(three_dimensional)
    gmsh.options.set("nodes",std::string("xyz"));
  gmsh.write(*m,"rll.msh");

}

BOOST_AUTO_TEST_CASE( test_meshgen_ghost_at_end )
{
  SharedPtr<grid::Grid> grid(grid::Grid::create("O8"));

  atlas::util::Config cfg;
  cfg.set("part",1);
  cfg.set("nb_parts",8);
  SharedPtr<MeshGenerator> meshgenerator( new ReducedGridMeshGenerator(cfg) );
  SharedPtr<mesh::Mesh> mesh ( meshgenerator->generate(*grid) );
  const array::ArrayView<int,1> part( mesh->nodes().partition() );
  const array::ArrayView<int,1> ghost( mesh->nodes().ghost() );
  const array::ArrayView<int,1> flags( mesh->nodes().field("flags") );

  Log::info() << "partition = [ ";
  for( size_t jnode=0; jnode<part.size(); ++jnode )
  {
    Log::info() << part(jnode) << " ";
  }
  Log::info() << "]" << std::endl;

  Log::info() << "ghost     = [ ";
  for( size_t jnode=0; jnode<part.size(); ++jnode )
  {
    Log::info() << ghost(jnode) << " ";
  }
  Log::info() << "]" << std::endl;

  Log::info() << "flags     = [ ";
  for( size_t jnode=0; jnode<part.size(); ++jnode )
  {
    Log::info() << internals::Topology::check(flags(jnode),internals::Topology::GHOST) << " ";
    BOOST_CHECK_EQUAL( internals::Topology::check(flags(jnode),internals::Topology::GHOST), ghost(jnode) );
  }
  Log::info() << "]" << std::endl;

}

} // namespace test
} // namespace atlas
