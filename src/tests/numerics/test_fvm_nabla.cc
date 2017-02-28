/*
 * (C) Copyright 1996-2016 ECMWF.
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

#include "atlas/atlas.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/grid/Grid.h"
#include "atlas/internals/Parameters.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/numerics/Nabla.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "eckit/memory/ScopedPtr.h"
#include "eckit/memory/SharedPtr.h"


using namespace eckit;
using namespace atlas::numerics;
using namespace atlas::mesh::generators;
using namespace atlas::internals;
using namespace atlas::grid;

namespace atlas {
namespace test {

double dual_volume(const mesh::Mesh& mesh)
{
  const mesh::Nodes& nodes = mesh.nodes();
  int nb_nodes = nodes.size();
  const array::ArrayView<double,1> dual_volumes ( nodes.field("dual_volumes") );
  const array::ArrayView<gidx_t,1> glb_idx ( nodes.global_index() );
  const array::ArrayView<int,1>    is_ghost ( nodes.ghost() );
  double area=0;
  for( int node=0; node<nb_nodes; ++node )
  {
    if( ! is_ghost(node) )
    {
      area += dual_volumes(node);
    }
  }
  ECKIT_MPI_CHECK_RESULT( MPI_Allreduce( MPI_IN_PLACE, &area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ) );
  return area;
}


/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
void rotated_flow(const fvm::Method& fvm, field::Field& field, const double& beta)
{
  const double radius = fvm.radius();
  const double USCAL = 20.;
  const double pvel = USCAL/radius;
  const double deg2rad = M_PI/180.;

  array::ArrayView<double,2> lonlat_deg (fvm.mesh().nodes().lonlat());
  array::ArrayView<double,3> var (field);

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
void rotated_flow_magnitude(const fvm::Method& fvm, field::Field& field, const double& beta)
{
  const double radius = fvm.radius();
  const double USCAL = 20.;
  const double pvel = USCAL/radius;
  const double deg2rad = M_PI/180.;

  array::ArrayView<double,2> lonlat_deg (fvm.mesh().nodes().lonlat());
  array::ArrayView<double,2> var (field);

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

struct AtlasFixture {
    AtlasFixture()
    {
      atlas_init(boost::unit_test::framework::master_test_suite().argc,
                 boost::unit_test::framework::master_test_suite().argv);
    }
    ~AtlasFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_factory )
{
  BOOST_CHECK( NablaFactory::has("fvm") );
}

BOOST_AUTO_TEST_CASE( test_build )
{
  Log::info() << "test_build" << std::endl;
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("Structured") );
  SharedPtr<mesh::Mesh> mesh( meshgenerator->generate( grid::Grid("O16") ) );
  const double R = util::Earth::radiusInMeters();
  fvm::Method fvm(*mesh,util::Config("radius",R));
  SharedPtr<Nabla> nabla ( Nabla::create(fvm) );

  double spherical_area = 360.*180.;
  BOOST_CHECK_CLOSE(dual_volume(*mesh),spherical_area,5.0);

}


BOOST_AUTO_TEST_CASE( test_grad )
{
  Log::info() << "test_grad" << std::endl;
  const double radius = util::Earth::radiusInMeters();
//  const double radius = 1.;
  grid::Grid grid(griduid());
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("Structured") );
  SharedPtr<mesh::Mesh> mesh( meshgenerator->generate(grid) );
  fvm::Method fvm(*mesh, util::Config("radius",radius));
  SharedPtr<Nabla> nabla ( Nabla::create(fvm) );

  array::ArrayView<double,2> lonlat( mesh->nodes().lonlat() );
  size_t nnodes = mesh->nodes().size();
  size_t nlev = 1;

  field::FieldSet fields;
  fields.add( fvm.node_columns().createField<double>("scalar",nlev,array::make_shape(1)) );
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

  nabla->gradient(fields["scalar"],fields["grad"]);
  nabla->gradient(fields["rscalar"],fields["rgrad"]);
  array::ArrayView<double,2> xder( fields["xder"] );
  array::ArrayView<double,2> yder( fields["yder"] );
  array::ArrayView<double,2> rxder( fields["rxder"] );
  array::ArrayView<double,2> ryder( fields["ryder"] );
  const array::ArrayView<double,3> grad( fields["grad"] );
  const array::ArrayView<double,3> rgrad( fields["rgrad"] );
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
    output::Gmsh(grid.name()+".msh").write(*mesh);
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
  grid::Grid grid(griduid());
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("Structured") );
  SharedPtr<mesh::Mesh> mesh( meshgenerator->generate(grid) );
  fvm::Method fvm(*mesh, util::Config("radius",radius));
  SharedPtr<Nabla> nabla ( Nabla::create(fvm) );

  array::ArrayView<double,2> lonlat( mesh->nodes().lonlat() );
  size_t nlev = 1;

  field::FieldSet fields;
  fields.add( fvm.node_columns().createField<double>("wind",nlev,array::make_shape(2)) );
  fields.add( fvm.node_columns().createField<double>("div",nlev) );

  rotated_flow(fvm,fields["wind"],M_PI_2*0.75);

  nabla->divergence(fields["wind"],fields["div"]);

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
  grid::Grid grid(griduid());
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("Structured") );
  SharedPtr<mesh::Mesh> mesh( meshgenerator->generate(grid) );
  fvm::Method fvm(*mesh, util::Config("radius",radius));
  SharedPtr<Nabla> nabla ( Nabla::create(fvm) );

  array::ArrayView<double,2> lonlat( mesh->nodes().lonlat() );
  size_t nlev = 1;

  field::FieldSet fields;
  fields.add( fvm.node_columns().createField<double>("wind",nlev,array::make_shape(2)) );
  fields.add( fvm.node_columns().createField<double>("vor",nlev) );

  rotated_flow(fvm,fields["wind"],M_PI_2*0.75);

  nabla->curl(fields["wind"],fields["vor"]);

  fields.add( fvm.node_columns().createField<double>("windgrad",nlev,array::make_shape(2,2)));
  nabla->gradient(fields["wind"],fields["windgrad"]);

  fields.add( fvm.node_columns().createField<double>("windX") );
  fields.add( fvm.node_columns().createField<double>("windY") );
  fields.add( fvm.node_columns().createField<double>("windXgradX") );
  fields.add( fvm.node_columns().createField<double>("windXgradY") );
  fields.add( fvm.node_columns().createField<double>("windYgradX") );
  fields.add( fvm.node_columns().createField<double>("windYgradY") );
  array::ArrayView<double,3> wind(fields["wind"]);
  array::ArrayView<double,4> windgrad(fields["windgrad"]);

  array::ArrayView<double,1> windX(fields["windX"]);
  array::ArrayView<double,1> windY(fields["windY"]);
  array::ArrayView<double,1> windXgradX(fields["windXgradX"]);
  array::ArrayView<double,1> windXgradY(fields["windXgradY"]);
  array::ArrayView<double,1> windYgradX(fields["windYgradX"]);
  array::ArrayView<double,1> windYgradY(fields["windYgradY"]);
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
  grid::Grid grid(griduid());
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("Structured") );
  SharedPtr<mesh::Mesh> mesh( meshgenerator->generate(grid) );
  fvm::Method fvm(*mesh, util::Config("radius",radius));
  SharedPtr<Nabla> nabla ( Nabla::create(fvm) );

  array::ArrayView<double,2> lonlat( mesh->nodes().lonlat() );
  size_t nlev = 1;

  field::FieldSet fields;
  fields.add( fvm.node_columns().createField<double>("scal",nlev) );
  fields.add( fvm.node_columns().createField<double>("lapl",nlev) );

  rotated_flow_magnitude(fvm,fields["scal"],M_PI_2*0.75);

  nabla->laplacian(fields["scal"],fields["lapl"]);

  // output to gmsh
  {
    fvm.node_columns().haloExchange(fields);
    output::Gmsh gmsh(grid.name()+"_fields.msh","a");
    gmsh.write(fields["lapl"]);
  }
}



} // namespace test
} // namespace atlas
