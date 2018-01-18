/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "eckit/types/Types.h"
#include "eckit/memory/ScopedPtr.h"
#include "atlas/library/Library.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/field/Field.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/output/Gmsh.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/MicroDeg.h"

#include "tests/AtlasTestEnvironment.h"
#include "eckit/testing/Test.h"


using namespace eckit;
using namespace eckit::testing;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_functionspace_StructuredColumns_no_halo" )
{
  int root=0;
  Grid grid("O8");
  util::Config config;
  config.set("halo",0);
  functionspace::StructuredColumns fs(grid, grid::Partitioner("equal_regions"), config );

  Field field     = fs.createField<double>( option::name("field") );
  Field field_glb = fs.createField<double>( option::name("field_global") | option::global(root));

  auto value     = array::make_view<double,1>( field );
  auto value_glb = array::make_view<double,1>( field_glb );

  value.assign(parallel::mpi::comm().rank());

  fs.gather(field,field_glb);

  Log::info() << "field checksum = " << fs.checksum(field) << std::endl;

//  for( size_t j=0; j<value_glb.size(); ++j )
//    Log::info() << value_glb(j) << " ";
//  Log::info() << std::endl;

  if( parallel::mpi::comm().rank() == root && parallel::mpi::comm().size() == 5 ) {
    std::vector<double> check{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

    EXPECT( value_glb.size() == check.size() );
    for( size_t j=0; j<value_glb.size(); ++j ) {
      EXPECT( value_glb(j) == check[j] );
    }
  }

  output::Gmsh gmsh("structured.msh");

  gmsh.write( MeshGenerator("structured").generate(grid) );
  gmsh.write( field );
}

CASE( "test_functionspace_StructuredColumns_halo" )
{
  ATLAS_DEBUG_VAR( parallel::mpi::comm().size() );
  int root=0;
//  grid::StructuredGrid grid(
//      grid::StructuredGrid::XSpace( {0.,360.} , {2,4,6,6,4,2} , false ),
//      grid::StructuredGrid::YSpace( grid::LinearSpacing( {80.,-80.}, 6 ) ),
//      Projection(),
//      Domain() );

  std::string gridname = eckit::Resource<std::string>("--grid","O8");

  grid::StructuredGrid grid(gridname);

  int halo = eckit::Resource<int>("--halo",2);
  util::Config config;
  config.set("halo",halo);
  functionspace::StructuredColumns fs(grid, grid::Partitioner("equal_regions"), config );

  Field field     = fs.createField<long>( option::name("field") );

  auto value = array::make_view<long,1>( field );
  auto xy    = array::make_view<double,2>( fs.xy() );
  auto g     = array::make_view<gidx_t,1>( fs.global_index() );
  auto r     = array::make_view<idx_t,1>( fs.remote_index() );
  auto p     = array::make_view<int,1>( fs.partition() );

  for( idx_t j=fs.j_begin(); j<fs.j_end(); ++j ) {
    for( idx_t i=fs.i_begin(j); i<fs.i_end(j); ++i ) {
      idx_t n = fs.index(i,j);
      value(n) = util::microdeg(xy(n,XX));
    }
  }

  //EXPECT( fs.checksum(field) == "cef2694016492d408fa157b7c59ce741" );

  fs.haloExchange(field);

  //EXPECT( fs.checksum(field) == "cef2694016492d408fa157b7c59ce741" );

  eckit::PathName filepath("test_functionspace_StructuredColumns_halo_p"+std::to_string(parallel::mpi::comm().rank())+".py");

  std::ofstream f(filepath.asString().c_str(), std::ios::trunc );

  f << "\n" "import matplotlib.pyplot as plt"
       "\n" "from matplotlib.path import Path"
       "\n" "import matplotlib.patches as patches"
       "\n" ""
       "\n" "from itertools import cycle"
       "\n" "import matplotlib.cm as cm"
       "\n" "import numpy as np"
       "\n" ""
       "\n" "fig = plt.figure(figsize=(20,10))"
       "\n" "ax = fig.add_subplot(111,aspect='equal')"
       "\n" "";

  double xmin= std::numeric_limits<double>::max();
  double xmax=-std::numeric_limits<double>::max();
  double ymin= std::numeric_limits<double>::max();
  double ymax=-std::numeric_limits<double>::max();
  f << "\n" "x = [";
  for( idx_t j=fs.j_begin_halo(); j<fs.j_end_halo(); ++j ) {
    for( idx_t i=fs.i_begin_halo(j); i<fs.i_end_halo(j); ++i ) {
      idx_t n = fs.index(i,j);
      f << xy(n,XX) << ", ";
      xmin = std::min(xmin,xy(n,XX));
      xmax = std::max(xmax,xy(n,XX));
    }
  }
  f << "]";

  f << "\n" "y = [";
  for( idx_t j=fs.j_begin_halo(); j<fs.j_end_halo(); ++j ) {
    for( idx_t i=fs.i_begin_halo(j); i<fs.i_end_halo(j); ++i ) {
      idx_t n = fs.index(i,j);
      f << xy(n,YY) << ", ";
      ymin = std::min(ymin,xy(n,YY));
      ymax = std::max(ymax,xy(n,YY));
    }
  }
  f << "]";

  f << "\n" "g = [";
  for( idx_t j=fs.j_begin_halo(); j<fs.j_end_halo(); ++j ) {
    for( idx_t i=fs.i_begin_halo(j); i<fs.i_end_halo(j); ++i ) {
      idx_t n = fs.index(i,j);
      f << g(n) << ", ";
    }
  }
  f << "]";

  f << "\n" "p = [";
  for( idx_t j=fs.j_begin_halo(); j<fs.j_end_halo(); ++j ) {
    for( idx_t i=fs.i_begin_halo(j); i<fs.i_end_halo(j); ++i ) {
      idx_t n = fs.index(i,j);
      f << p(n) << ", ";
    }
  }
  f << "]";

  f << "\n" "r = [";
  for( idx_t j=fs.j_begin_halo(); j<fs.j_end_halo(); ++j ) {
    for( idx_t i=fs.i_begin_halo(j); i<fs.i_end_halo(j); ++i ) {
      idx_t n = fs.index(i,j);
      f << r(n) << ", ";
    }
  }
  f << "]";

  f << "\n" ""
       "\n" "c = [ cm.Paired( float(pp%13)/12. ) for pp in p ]"
       "\n" "ax.scatter(x, y, color=c, marker='o')"
       "\n" "for i in range("<<fs.size()<<"):"
       "\n" "  ax.annotate(g[i], (x[i],y[i]), fontsize=8)"
       "\n" "";
  f <<
       "\n" "ax.set_xlim( "<<std::min(0.,xmin)<<"-5, "<<std::max(360.,xmax)<<"+5)"
       "\n" "ax.set_ylim( "<<std::min(-90.,ymin)<<"-5, "<<std::max(90.,ymax)<<"+5)"
       "\n" "ax.set_xticks([0,45,90,135,180,225,270,315,360])"
       "\n" "ax.set_yticks([-90,-45,0,45,90])"
       "\n" "plt.grid()"
       "\n" "plt.show()"
       "\n";
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char **argv) {
    atlas::test::AtlasTestEnvironment env( argc, argv );
    return run_tests ( argc, argv, false );
}
