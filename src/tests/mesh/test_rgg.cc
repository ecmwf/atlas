/*
 * (C) Copyright 1996-2017 ECMWF.
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

#include "atlas/atlas.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/internals/atlas_config.h"
#include "atlas/grid/detail/spacing/gaussian/Latitudes.h"
#include "atlas/grid.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/grid/detail/partitioners/EqualRegionsPartitioner.h"
#include "atlas/grid/detail/grid/reduced/ReducedGaussian.h"
#include "atlas/grid/detail/grid/CustomStructured.h"
#include "atlas/output/Gmsh.h"
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
#include "atlas/internals/Parameters.h"
#include "atlas/util/Config.h"
#include "atlas/internals/Bitflags.h"
#include "atlas/internals/Debug.h"

#include "tests/AtlasFixture.h"


namespace atlas {
namespace grid {
namespace spacing {
namespace gaussian {
void compute_gaussian_quadrature_npole_equator(const size_t N, double lats[], double weights[]);
}
}
}
}

#define DISABLE if(0)
#define ENABLE if(1)

namespace atlas {
namespace test {

using eckit::geometry::LAT;
using eckit::geometry::LON;

class DebugGrid: public grid::detail::grid::reduced::ReducedGaussian { public: DebugGrid(); };
DebugGrid::DebugGrid()
{
  int N=5;
  long lon[] = {
    6,
    10,
    18,
    22,
    22,
  };
  grid::detail::grid::reduced::ReducedGaussian::setup(N,lon);
}

static grid::StructuredGrid debug_grid() { return grid::StructuredGrid( new DebugGrid() ); }


class MinimalGrid:   public grid::detail::grid::reduced::ReducedGaussian {
	public:
		MinimalGrid(int N, long lon[]) : ReducedGaussian(N,lon) {}
};

static grid::StructuredGrid minimal_grid(int N, long lon[]) { return grid::StructuredGrid( new MinimalGrid(N,lon) ); }


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

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_eq_caps )
{
  std::vector<int>    n_regions;
  std::vector<double> s_cap;

  grid::detail::partitioners::eq_caps(6, n_regions, s_cap);
  BOOST_CHECK_EQUAL( n_regions.size(), 3 );
  BOOST_CHECK_EQUAL( n_regions[0], 1 );
  BOOST_CHECK_EQUAL( n_regions[1], 4 );
  BOOST_CHECK_EQUAL( n_regions[2], 1 );

  grid::detail::partitioners::eq_caps(10, n_regions, s_cap);
  BOOST_CHECK_EQUAL( n_regions.size(), 4 );
  BOOST_CHECK_EQUAL( n_regions[0], 1 );
  BOOST_CHECK_EQUAL( n_regions[1], 4 );
  BOOST_CHECK_EQUAL( n_regions[2], 4 );
  BOOST_CHECK_EQUAL( n_regions[3], 1 );
}

BOOST_AUTO_TEST_CASE( test_partitioner )
{
  grid::Grid g( "S4x2" );

  // 12 partitions
  {
    grid::detail::partitioners::EqualRegionsPartitioner partitioner(g,12);
    BOOST_CHECK_EQUAL( partitioner.nb_bands(),    4 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(0), 1 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(1), 5 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(2), 5 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(3), 1 );
  }

  // 24 partitions
  {
    grid::detail::partitioners::EqualRegionsPartitioner partitioner(g,24);
    BOOST_CHECK_EQUAL( partitioner.nb_bands(),     5 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(0),  1 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(1),  6 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(2), 10 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(3),  6 );
    BOOST_CHECK_EQUAL( partitioner.nb_regions(4),  1 );
  }

  // 48 partitions
  {
    grid::detail::partitioners::EqualRegionsPartitioner partitioner(g,48);
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
    grid::detail::partitioners::EqualRegionsPartitioner partitioner(g,96);
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

BOOST_AUTO_TEST_CASE( test_gaussian_latitudes )
{
  std::vector< double > factory_latitudes;
  std::vector< double > computed_latitudes;
  std::vector< double > computed_weights;


  size_t size_test_N = 19;

  size_t test_N[] = {16,24,32,48,64,80,96,128,160,
                     200,256,320,400,512,576,640,
                     800,1024,1280,1600,2000,4000,8000};

  for( size_t i=0; i<size_test_N; ++i )
  {
    size_t N = test_N[i];
    Log::info() << "Testing gaussian latitude " << N << std::endl;
    factory_latitudes.resize(N);
    computed_latitudes.resize(N);
    computed_weights.resize(N);
    //grid::gaussian::latitudes::gaussian_latitudes_npole_equator (N, factory_latitudes.data());
    //grid::gaussian::latitudes::compute_gaussian_quadrature_npole_equator(N, computed_latitudes.data(), computed_weights.data());
    grid::spacing::gaussian::gaussian_latitudes_npole_equator (N, factory_latitudes.data());
    grid::spacing::gaussian::compute_gaussian_quadrature_npole_equator(N, computed_latitudes.data(), computed_weights.data());
    double wsum=0;
    for( size_t i=0; i<N; ++i )
    {
      BOOST_CHECK_CLOSE( computed_latitudes[i] , factory_latitudes[i], 0.0000001 );
      wsum += computed_weights[i];
    }
    BOOST_CHECK_CLOSE( wsum*2. , 1. , 0.0000001 );
  }
}

BOOST_AUTO_TEST_CASE( test_rgg_meshgen_one_part )
{
  mesh::Mesh* m;
  util::Config default_opts;
  default_opts.set("nb_parts",1);
  default_opts.set("part",0);
//  generate.options.set("nb_parts",1);
//  generate.options.set("part",    0);
DISABLE{  // This is all valid for meshes generated with MINIMAL NB TRIAGS
  ENABLE {

    mesh::generators::Structured generate(
          default_opts
          ("3d",true)
          ("include_pole",false)
          );
    m = generate( atlas::test::debug_grid() );
    BOOST_CHECK_EQUAL( m->nodes().size(), 156 );
    BOOST_CHECK_EQUAL( m->cells().elements(0).size(), 134 );
    BOOST_CHECK_EQUAL( m->cells().elements(1).size(),  32 );
//    BOOST_CHECK_EQUAL( m->nodes().metadata().get<size_t>("nb_owned"),    156 );
//    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("nb_owned"),    134 );
//    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("nb_owned"),     32 );
    delete m;
  }

  ENABLE {
    mesh::generators::Structured generate(
          default_opts
          ("3d",false)
          ("include_pole",false)
          );
    m = generate( atlas::test::debug_grid() );
    BOOST_CHECK_EQUAL( m->nodes().size(), 166 );
    BOOST_CHECK_EQUAL( m->cells().elements(0).size(), 134 );
    BOOST_CHECK_EQUAL( m->cells().elements(1).size(),  32 );
//    BOOST_CHECK_EQUAL( m->nodes().metadata().get<size_t>("nb_owned"),    166 );
//    BOOST_CHECK_EQUAL( m->function_space("quads" ).metadata().get<size_t>("nb_owned"),    134 );
//    BOOST_CHECK_EQUAL( m->function_space("triags").metadata().get<size_t>("nb_owned"),     32 );
    delete m;
  }

  ENABLE {
    mesh::generators::Structured generate(
          default_opts
          ("3d",true)
          ("include_pole",true)
          );
    m = generate( atlas::test::debug_grid() );
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
    mesh::generators::Structured generate(
          default_opts
          ("3d",false)
          ("include_pole",false)
          );
    int nlat=2;
    long lon[] = { 4, 6 };
    mesh = generate( test::minimal_grid(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->nodes().size(), 24 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(0).size(), 14 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(1).size(),  4 );

    double max_lat = test::minimal_grid(nlat,lon).y().front();
    BOOST_CHECK_CLOSE( test::compute_lonlat_area(*mesh), 2.*M_PI*2.*max_lat, 1e-8 );
    output::Gmsh("minimal2.msh").write(*mesh);
    delete mesh;
  }
  // 3 latitudes
  ENABLE {
    mesh::generators::Structured generate(
          default_opts
          ("3d",false)
          ("include_pole",false)
          );
    int nlat=3;
    long lon[] = { 4, 6, 8 };
    mesh = generate( test::minimal_grid(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->nodes().size(), 42 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(0).size(), 28 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(1).size(),  8 );
    output::Gmsh("minimal3.msh").write(*mesh);
    delete mesh;
  }
  // 4 latitudes
  ENABLE {
    mesh::generators::Structured generate(
          default_opts
          ("3d",false)
          ("include_pole",false)
          );
    int nlat=4;
    long lon[] = { 4, 6, 8, 10 };
    mesh = generate( test::minimal_grid(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->nodes().size(), 64 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(0).size(), 46 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(1).size(), 12 );
    output::Gmsh("minimal4.msh").write(*mesh);
    delete mesh;
  }
  // 5 latitudes WIP
  ENABLE {
    mesh::generators::Structured generate(
          default_opts
          ("3d",false)
          ("include_pole",false)
          );
    int nlat=5;
    long lon[] = { 6, 10, 18, 22, 22 };
    mesh = generate( test::minimal_grid(nlat,lon) );
    BOOST_CHECK_EQUAL( mesh->nodes().size(), 166 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(0).size(), 134 );
    BOOST_CHECK_EQUAL( mesh->cells().elements(1).size(),  32 );
    output::Gmsh("minimal5.msh").write(*mesh);
    delete mesh;
  }
}
}

BOOST_AUTO_TEST_CASE( test_rgg_meshgen_many_parts )
{

  BOOST_CHECK( grid::detail::partitioners::PartitionerFactory::has("EqualRegions") );
  size_t nb_parts = 20;
          //  Alternative grid for debugging
          //  int nlat=10;
          //  long lon[] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 };
          //  test::MinimalGrid grid(nlat,lon);
  grid::StructuredGrid grid = grid::Grid("N32");
  //RegularGrid grid(128,64);

  /*
  std::cout << grid.spec() << std::endl;
  for (int jlat=0;jlat<2*nlat; jlat++) {
  	std::cout << grid.lon(jlat,0) << ", ... , " << "," << grid.lon(jlat,9) << grid.lon(jlat,10) << std::endl;
  }
  ASSERT(0);
  */

  double max_lat = grid.y().front();
  double check_area = 360.*2.*max_lat;
  double area = 0;
  int nodes[]  = {313,332,336,338,334,337,348,359,360,361,360,360,359,370,321,334,338,335,348,315};
  int quads[]  = {243,290,293,294,291,293,310,320,321,321,320,321,320,331,278,291,294,293,305,245};
  int triags[] = { 42, 13, 12, 13, 12, 14,  0,  1,  0,  1,  1,  0,  1,  0, 14, 12, 13, 11, 14, 42};
  int nb_owned = 0;

  std::vector<int> all_owned    ( grid.npts()+grid.ny()+1, -1 );

  for(size_t p = 0; p < nb_parts; ++p)
  {

    DEBUG_VAR(p);

    mesh::generators::Structured generate ( util::Config
           ("partitioner","EqualRegions")
           ("nb_parts",nb_parts)
           ("part",p)
           ("include_pole",false)
           ("3d",false) );

    mesh::Mesh::Ptr m( generate( grid ) );
    DEBUG_HERE();
    m->metadata().set("part",p);
    BOOST_TEST_CHECKPOINT("generated grid " << p);
    array::ArrayView<int,1> part( m->nodes().partition() );
    array::ArrayView<gidx_t,1> gidx( m->nodes().global_index() );

    area += test::compute_lonlat_area(*m);
    DEBUG_HERE();

    DISABLE {  // This is all valid for meshes generated with MINIMAL NB TRIAGS
    if( nb_parts == 20 )
    {
      BOOST_CHECK_EQUAL( m->nodes().size(), nodes[p]  );
      BOOST_CHECK_EQUAL( m->cells().elements(0).size(), quads[p]  );
      BOOST_CHECK_EQUAL( m->cells().elements(1).size(), triags[p] );
    }
    }
    DEBUG_HERE();

    output::Gmsh("T63.msh").write(*m);


    mesh::Nodes& nodes = m->nodes();
    size_t nb_nodes = nodes.size();

    // Test if all nodes are connected
    {
      std::vector<int> node_elem_connections( nb_nodes, 0 );

      const mesh::HybridElements::Connectivity& cell_node_connectivity = m->cells().node_connectivity();
      for( size_t jelem=0; jelem<m->cells().size(); ++jelem )
      {
        for( size_t jnode=0; jnode<cell_node_connectivity.cols(jelem); ++jnode )
          node_elem_connections[ cell_node_connectivity(jelem,jnode) ]++;
      }
      for( size_t jnode=0; jnode<nb_nodes; ++jnode )
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
    for( size_t n=0; n<nb_nodes; ++n )
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
  BOOST_CHECK_EQUAL( nb_owned, grid.npts()+grid.ny() );

  BOOST_CHECK_CLOSE( area, check_area, 1e-10 );

}

BOOST_AUTO_TEST_CASE( test_reduced_lonlat )
{
  DEBUG_HERE();
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
  grid::StructuredGrid grid( new grid::detail::grid::CustomStructured(N,lat,lon) );

  bool three_dimensional = true;

  mesh::generators::Structured generate( util::Config
        ("3d",three_dimensional)
        ("triangulate",false) );

  mesh::Mesh::Ptr m (generate(grid));

  util::Config options;
  if( three_dimensional ) options.set("coordinates","xyz");
  output::Gmsh gmsh("rll.msh",options);
  gmsh.write(*m);

}

BOOST_AUTO_TEST_CASE( test_meshgen_ghost_at_end )
{
  DEBUG_HERE();

  grid::Grid grid("O8");

  atlas::util::Config cfg;
  cfg.set("part",1);
  cfg.set("nb_parts",8);
  eckit::SharedPtr<mesh::generators::MeshGenerator> meshgenerator( new mesh::generators::Structured(cfg) );
  eckit::SharedPtr<mesh::Mesh> mesh ( meshgenerator->generate(grid) );
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
