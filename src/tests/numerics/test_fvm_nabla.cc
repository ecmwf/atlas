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
#include "eckit/memory/ScopedPtr.h"
#include "eckit/memory/SharedPtr.h"
#include "atlas/atlas.h"
#include "atlas/numerics/Nabla.h"
#include "atlas/util/Config.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/numerics/fvm/Method.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/field/Field.h"
#include "atlas/internals/Parameters.h"
#include "atlas/util/Constants.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/field/FieldSet.h"
#include "atlas/parallel/mpi/mpi.h"

using namespace eckit;
using namespace atlas::numerics;
using namespace atlas::mesh::generators;
using namespace atlas::internals;
using namespace atlas::grid;

namespace atlas {
namespace test {

double dual_volume(mesh::Mesh& mesh)
{
  mesh::Nodes& nodes = mesh.nodes();
  int nb_nodes = nodes.size();
  const array::ArrayView<double,1> dual_volumes ( nodes.field("dual_volumes") );
  const array::ArrayView<gidx_t,1> glb_idx ( nodes.global_index() );
  const array::ArrayView<int,1> is_ghost ( nodes.ghost() );
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



struct AtlasFixture {
    AtlasFixture()  { atlas_init(); }
    ~AtlasFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_factory )
{
  BOOST_CHECK( NablaFactory::has("fvm") );
}

BOOST_AUTO_TEST_CASE( test_build )
{
  SharedPtr<Grid> grid ( Grid::create("O32") );
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("Structured") );
  SharedPtr<mesh::Mesh> mesh( meshgenerator->generate(*grid) );
  const double R = util::Earth::radiusInMeters();
  fvm::Method fvm(*mesh,util::Config("radius",R));
  SharedPtr<Nabla> nabla ( Nabla::create(fvm) );

  double spherical_area = 4.*M_PI*R*R;
  BOOST_CHECK_CLOSE(dual_volume(*mesh),spherical_area,0.001*spherical_area);

}


BOOST_AUTO_TEST_CASE( test_grad )
{
  const double radius = util::Earth::radiusInMeters();
//  const double radius = 1.;
  SharedPtr<Grid> grid ( Grid::create("O32") );
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("Structured") );
  SharedPtr<mesh::Mesh> mesh( meshgenerator->generate(*grid) );
  fvm::Method fvm(*mesh, util::Config("radius",radius));
  SharedPtr<Nabla> nabla ( Nabla::create(fvm) );

  array::ArrayView<double,2> lonlat( mesh->nodes().lonlat() );
  size_t nnodes = mesh->nodes().size();
  size_t nlev = 1;

  field::FieldSet fields;
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
    util::io::Gmsh().write(*mesh,grid->shortName()+".msh");
    util::io::Gmsh().write(fields["scalar"],grid->shortName()+"_fields.msh");
    util::io::Gmsh().write(fields["xder"],grid->shortName()+"_fields.msh",std::ios::app);
    util::io::Gmsh().write(fields["yder"],grid->shortName()+"_fields.msh",std::ios::app);
    util::io::Gmsh().write(fields["rscalar"],grid->shortName()+"_fields.msh",std::ios::app);
    util::io::Gmsh().write(fields["rxder"],grid->shortName()+"_fields.msh",std::ios::app);
    util::io::Gmsh().write(fields["ryder"],grid->shortName()+"_fields.msh",std::ios::app);

    //    util::io::Gmsh().write(fields["exact_yder"],grid->shortName()+"_fields.msh",std::ios::app);
  }
}


BOOST_AUTO_TEST_CASE( test_div )
{
  const double radius = util::Earth::radiusInMeters();
//  const double radius = 1.;
  SharedPtr<Grid> grid ( Grid::create("O32") );
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("Structured") );
  SharedPtr<mesh::Mesh> mesh( meshgenerator->generate(*grid) );
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
    util::io::Gmsh().write(*mesh,grid->shortName()+".msh");
    util::io::Gmsh().write(fields["wind"],grid->shortName()+"_fields.msh",std::ios::app);
    util::io::Gmsh().write(fields["div"],grid->shortName()+"_fields.msh",std::ios::app);
  }
}

BOOST_AUTO_TEST_CASE( test_curl )
{
  const double radius = util::Earth::radiusInMeters();
//  const double radius = 1.;
  SharedPtr<Grid> grid ( Grid::create("O32") );
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("Structured") );
  SharedPtr<mesh::Mesh> mesh( meshgenerator->generate(*grid) );
  fvm::Method fvm(*mesh, util::Config("radius",radius));
  SharedPtr<Nabla> nabla ( Nabla::create(fvm) );

  array::ArrayView<double,2> lonlat( mesh->nodes().lonlat() );
  size_t nlev = 1;

  field::FieldSet fields;
  fields.add( fvm.node_columns().createField<double>("wind",nlev,array::make_shape(2)) );
  fields.add( fvm.node_columns().createField<double>("vor",nlev) );

  rotated_flow(fvm,fields["wind"],M_PI_2*0.75);

  nabla->curl(fields["wind"],fields["vor"]);

  // output to gmsh
  {
    fvm.node_columns().haloExchange(fields);
    util::io::Gmsh().write(*mesh,grid->shortName()+".msh");
//    util::io::Gmsh().write(fields["wind"],grid->shortName()+"_fields.msh",std::ios::app);
    util::io::Gmsh().write(fields["vor"],grid->shortName()+"_fields.msh",std::ios::app);
  }
}

BOOST_AUTO_TEST_CASE( test_lapl )
{
  const double radius = util::Earth::radiusInMeters();
//  const double radius = 1.;
  SharedPtr<Grid> grid ( Grid::create("O32") );
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("Structured") );
  SharedPtr<mesh::Mesh> mesh( meshgenerator->generate(*grid) );
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
    util::io::Gmsh().write(*mesh,grid->shortName()+".msh");
//    util::io::Gmsh().write(fields["wind"],grid->shortName()+"_fields.msh",std::ios::app);
    util::io::Gmsh().write(fields["lapl"],grid->shortName()+"_fields.msh",std::ios::app);
  }
}



} // namespace test
} // namespace atlas
