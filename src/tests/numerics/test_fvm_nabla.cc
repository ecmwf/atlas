/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestNablaEdgeBasedFiniteVolume
#include "ecbuild/boost_test_framework.h"

#include <cmath>
#include <iostream>

#include "atlas/library/Library.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Grid.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/numerics/Nabla.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Earth.h"
#include "atlas/grid/Partitioner.h"

#include "tests/AtlasFixture.h"


using namespace eckit;
using namespace atlas::numerics;
using namespace atlas::meshgenerator;
using namespace atlas::grid;

namespace atlas {
namespace test {

double dual_volume(const Mesh& mesh)
{
  const mesh::Nodes& nodes = mesh.nodes();
  int nb_nodes = nodes.size();
  const array::ArrayView<double,1> dual_volumes = array::make_view<double,1>( nodes.field("dual_volumes") );
  const array::ArrayView<int   ,1> is_ghost     = array::make_view<int   ,1>( nodes.ghost() );
  double area=0;
  for( int node=0; node<nb_nodes; ++node )
  {
    if( ! is_ghost(node) )
    {
      area += dual_volumes(node);
    }
  }

  parallel::mpi::comm().allReduceInPlace(area, eckit::mpi::sum());

  return area;
}


/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
void rotated_flow(const fvm::Method& fvm, Field& field, const double& beta)
{
  const double radius = fvm.radius();
  const double USCAL = 20.;
  const double pvel = USCAL/radius;
  const double deg2rad = M_PI/180.;

  array::ArrayView<double,2> lonlat_deg = array::make_view<double,2>(fvm.mesh().nodes().lonlat());
  array::ArrayView<double,3> var        = array::make_view<double,3>(field);

  size_t nnodes = fvm.mesh().nodes().size();
  for( size_t jnode=0; jnode<nnodes; ++jnode )
  {
     double x = lonlat_deg(jnode,LON) * deg2rad;
     double y = lonlat_deg(jnode,LAT) * deg2rad;
     double Ux =  pvel*(std::cos(beta)+std::tan(y)*std::cos(x)*std::sin(beta))*radius*std::cos(y);
     double Uy = -pvel*std::sin(x)*std::sin(beta)*radius;
     for( size_t jlev=0; jlev<field.levels(); ++jlev)
     {
       var(jnode,jlev,LON) = Ux;
       var(jnode,jlev,LAT) = Uy;
     }
  }
}

/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
void rotated_flow_magnitude(const fvm::Method& fvm, Field& field, const double& beta)
{
  const double radius = fvm.radius();
  const double USCAL = 20.;
  const double pvel = USCAL/radius;
  const double deg2rad = M_PI/180.;

  array::ArrayView<double,2> lonlat_deg = array::make_view<double,2> (fvm.mesh().nodes().lonlat());
  array::ArrayView<double,2> var        = array::make_view<double,2> (field);

  size_t nnodes = fvm.mesh().nodes().size();
  for( size_t jnode=0; jnode<nnodes; ++jnode )
  {
     double x = lonlat_deg(jnode,LON) * deg2rad;
     double y = lonlat_deg(jnode,LAT) * deg2rad;
     double Ux =  pvel*(std::cos(beta)+std::tan(y)*std::cos(x)*std::sin(beta))*radius*std::cos(y);
     double Uy = -pvel*std::sin(x)*std::sin(beta)*radius;
     for( size_t jlev=0; jlev<field.levels(); ++jlev )
       var(jnode,jlev) = std::sqrt(Ux*Ux+Uy*Uy);
  }
}

static std::string griduid() { return "Slat80"; }


BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_factory )
{
  BOOST_CHECK( NablaFactory::has("fvm") );
}

BOOST_AUTO_TEST_CASE( test_build )
{
  Log::info() << "test_build" << std::endl;
  MeshGenerator meshgenerator ("structured" );
  Mesh mesh = meshgenerator.generate( Grid("O16") );
  const double R = util::Earth::radiusInMeters();
  fvm::Method fvm(mesh,util::Config("radius",R));
  Nabla nabla( fvm );

  double spherical_area = 360.*180.;
  BOOST_CHECK_CLOSE(dual_volume(mesh),spherical_area,5.0);

}


BOOST_AUTO_TEST_CASE( test_grad )
{
  Log::info() << "test_grad" << std::endl;
  const double radius = util::Earth::radiusInMeters();
//  const double radius = 1.;
  Grid grid(griduid());
  MeshGenerator meshgenerator("structured");
  Mesh mesh = meshgenerator.generate(grid, Distribution(grid,Partitioner("equal_regions")) );
  fvm::Method fvm(mesh, util::Config("radius",radius));
  Nabla nabla( fvm );

  size_t nnodes = mesh.nodes().size();
  size_t nlev = 1;

  FieldSet fields;
  fields.add( fvm.node_columns().createField<double>("scalar",nlev) );
  fields.add( fvm.node_columns().createField<double>("rscalar",nlev) );
  fields.add( fvm.node_columns().createField<double>("grad",nlev,array::make_shape(2)) );
  fields.add( fvm.node_columns().createField<double>("rgrad",nlev,array::make_shape(2)) );
  fields.add( fvm.node_columns().createField<double>("xder",nlev) );
  fields.add( fvm.node_columns().createField<double>("yder",nlev) );
  fields.add( fvm.node_columns().createField<double>("rxder",nlev) );
  fields.add( fvm.node_columns().createField<double>("ryder",nlev) );

  //  fields.add( fvm.createField<double>("exact_yder",nlev) );

//  const double deg2rad = M_PI/180.;
//  array::ArrayView<double,2> var( fields["scalar"] );
////  array::ArrayView<double,2> exact_yder( fields["exact_yder"] );
//  for( size_t jnode=0; jnode< nnodes ; ++jnode )
//  {
//    const double y  = lonlat(jnode,LAT) * deg2rad ;

//    for(size_t jlev = 0; jlev < nlev; ++jlev) {
//      var(jnode,jlev)        = std::sin(4.*y);
////      exact_yder(jnode,jlev) = 4.*std::cos(4.*y)/radius;
//    }
//  }

  rotated_flow_magnitude(fvm,fields["scalar"],0.);
  rotated_flow_magnitude(fvm,fields["rscalar"],M_PI_2*0.75);

  nabla.gradient(fields["scalar"],fields["grad"]);
  nabla.gradient(fields["rscalar"],fields["rgrad"]);
  array::ArrayView<double,2> xder  = array::make_view<double,2>( fields["xder"] );
  array::ArrayView<double,2> yder  = array::make_view<double,2>( fields["yder"] );
  array::ArrayView<double,2> rxder = array::make_view<double,2>( fields["rxder"] );
  array::ArrayView<double,2> ryder = array::make_view<double,2>( fields["ryder"] );
  const array::ArrayView<double,3> grad  = array::make_view<double,3>( fields["grad"] );
  const array::ArrayView<double,3> rgrad = array::make_view<double,3>( fields["rgrad"] );
  for( size_t jnode=0; jnode< nnodes ; ++jnode )
  {
    for(size_t jlev = 0; jlev < nlev; ++jlev) {
      xder(jnode,jlev) = grad(jnode,jlev,LON);
      yder(jnode,jlev) = grad(jnode,jlev,LAT);
      rxder(jnode,jlev) = rgrad(jnode,jlev,LON);
      ryder(jnode,jlev) = rgrad(jnode,jlev,LAT);
    }
  }

  // output to gmsh
  {
    fvm.node_columns().haloExchange(fields);
    output::Gmsh(grid.name()+".msh").write(mesh);
    output::Gmsh gmsh_fields(grid.name()+"_fields.msh");
    gmsh_fields.write(fields["scalar"]);
    gmsh_fields.write(fields["xder"]);
    gmsh_fields.write(fields["yder"]);
    gmsh_fields.write(fields["rscalar"]);
    gmsh_fields.write(fields["rxder"]);
    gmsh_fields.write(fields["ryder"]);
  }
}


BOOST_AUTO_TEST_CASE( test_div )
{
  Log::info() << "test_div" << std::endl;
  const double radius = util::Earth::radiusInMeters();
//  const double radius = 1.;
  Grid grid(griduid());
  MeshGenerator meshgenerator("structured");
  Mesh mesh = meshgenerator.generate(grid, Distribution(grid,Partitioner("equal_regions")) );
  fvm::Method fvm(mesh, util::Config("radius",radius));
  Nabla nabla(fvm);

  size_t nlev = 1;

  FieldSet fields;
  fields.add( fvm.node_columns().createField<double>("wind",nlev,array::make_shape(2)) );
  fields.add( fvm.node_columns().createField<double>("div",nlev) );

  rotated_flow(fvm,fields["wind"],M_PI_2*0.75);

  nabla.divergence(fields["wind"],fields["div"]);

  // output to gmsh
  {
    fvm.node_columns().haloExchange(fields);
    output::Gmsh gmsh (grid.name()+"_fields.msh","a");
    gmsh.write(fields["wind"]);
    gmsh.write(fields["div"]);
  }
}

BOOST_AUTO_TEST_CASE( test_curl )
{
  Log::info() << "test_curl" << std::endl;
  const double radius = util::Earth::radiusInMeters();
//  const double radius = 1.;
  Grid grid(griduid());
  MeshGenerator meshgenerator("structured");
  Mesh mesh = meshgenerator.generate(grid, Distribution(grid,Partitioner("equal_regions")) );
  fvm::Method fvm(mesh, util::Config("radius",radius));
  Nabla nabla( fvm );

  size_t nlev = 1;

  FieldSet fields;
  fields.add( fvm.node_columns().createField<double>("wind",nlev,array::make_shape(2)) );
  fields.add( fvm.node_columns().createField<double>("vor",nlev) );

  rotated_flow(fvm,fields["wind"],M_PI_2*0.75);

  nabla.curl(fields["wind"],fields["vor"]);

  fields.add( fvm.node_columns().createField<double>("windgrad",nlev,array::make_shape(2,2)));
  nabla.gradient(fields["wind"],fields["windgrad"]);

  fields.add( fvm.node_columns().createField<double>("windX") );
  fields.add( fvm.node_columns().createField<double>("windY") );
  fields.add( fvm.node_columns().createField<double>("windXgradX") );
  fields.add( fvm.node_columns().createField<double>("windXgradY") );
  fields.add( fvm.node_columns().createField<double>("windYgradX") );
  fields.add( fvm.node_columns().createField<double>("windYgradY") );
  array::ArrayView<double,3> wind     = array::make_view<double,3>(fields["wind"]);
  array::ArrayView<double,4> windgrad = array::make_view<double,4>(fields["windgrad"]);

  array::ArrayView<double,1> windX      = array::make_view<double,1>(fields["windX"]);
  array::ArrayView<double,1> windY      = array::make_view<double,1>(fields["windY"]);
  array::ArrayView<double,1> windXgradX = array::make_view<double,1>(fields["windXgradX"]);
  array::ArrayView<double,1> windXgradY = array::make_view<double,1>(fields["windXgradY"]);
  array::ArrayView<double,1> windYgradX = array::make_view<double,1>(fields["windYgradX"]);
  array::ArrayView<double,1> windYgradY = array::make_view<double,1>(fields["windYgradY"]);
  for( size_t j=0; j<windX.size(); ++j )
  {
    windX(j) = wind(j,0,0);
    windY(j) = wind(j,0,1);
    windXgradX(j) = windgrad(j,0,0,0);
    windXgradY(j) = windgrad(j,0,0,1);
    windYgradX(j) = windgrad(j,0,1,0);
    windYgradY(j) = windgrad(j,0,1,1);
  }

  // output to gmsh
  {
    fvm.node_columns().haloExchange(fields);
    output::Gmsh gmsh(grid.name()+"_fields.msh","a");
    gmsh.write(fields["vor"]);
    gmsh.write(fields["windX"]);
    gmsh.write(fields["windXgradX"]);
    gmsh.write(fields["windXgradY"]);
    gmsh.write(fields["windY"]);
    gmsh.write(fields["windYgradX"]);
    gmsh.write(fields["windYgradY"]);
    gmsh.write(fields["windgrad"]);
  }





}

BOOST_AUTO_TEST_CASE( test_lapl )
{
  Log::info() << "test_lapl" << std::endl;
  const double radius = util::Earth::radiusInMeters();
//  const double radius = 1.;
  Grid grid(griduid());
  MeshGenerator meshgenerator("structured");
  Mesh mesh = meshgenerator.generate(grid, Distribution(grid,Partitioner("equal_regions")) );
  fvm::Method fvm(mesh, util::Config("radius",radius));
  Nabla nabla( fvm );

  size_t nlev = 1;

  FieldSet fields;
  fields.add( fvm.node_columns().createField<double>("scal",nlev) );
  fields.add( fvm.node_columns().createField<double>("lapl",nlev) );

  rotated_flow_magnitude(fvm,fields["scal"],M_PI_2*0.75);

  nabla.laplacian(fields["scal"],fields["lapl"]);

  // output to gmsh
  {
    fvm.node_columns().haloExchange(fields);
    output::Gmsh gmsh(grid.name()+"_fields.msh","a");
    gmsh.write(fields["lapl"]);
  }
}



} // namespace test
} // namespace atlas
