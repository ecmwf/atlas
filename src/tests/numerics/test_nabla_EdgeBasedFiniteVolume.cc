/*
 * (C) Copyright 1996-2015 ECMWF.
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
#include "atlas/Config.h"
#include "atlas/Grid.h"
#include "atlas/Mesh.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/functionspace/EdgeBasedFiniteVolume.h"
#include "atlas/Nodes.h"
#include "atlas/Field.h"
#include "atlas/Parameters.h"
#include "atlas/io/Gmsh.h"
#include "atlas/FieldSet.h"

using namespace eckit;
using namespace atlas::numerics;
using namespace atlas::meshgen;

namespace atlas {
namespace test {

struct AtlasFixture {
    AtlasFixture()  { atlas_init(); }
    ~AtlasFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( AtlasFixture )

BOOST_AUTO_TEST_CASE( test_factory )
{
  BOOST_CHECK( NablaFactory::has("EdgeBasedFiniteVolume") );
}

BOOST_AUTO_TEST_CASE( test_build )
{
  SharedPtr<Grid> grid ( Grid::create("oct.N24") );
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("ReducedGrid") );
  SharedPtr<Mesh> mesh( meshgenerator->generate(*grid) );
  functionspace::EdgeBasedFiniteVolume fvm(*mesh);
  SharedPtr<Nabla> nabla ( Nabla::create(fvm) );
}


BOOST_AUTO_TEST_CASE( test_grad )
{
  SharedPtr<Grid> grid ( Grid::create("oct.N24") );
  SharedPtr<MeshGenerator> meshgenerator ( MeshGenerator::create("ReducedGrid") );
  SharedPtr<Mesh> mesh( meshgenerator->generate(*grid) );
  functionspace::EdgeBasedFiniteVolume fvm(*mesh);
  SharedPtr<Nabla> nabla ( Nabla::create(fvm) );

  ArrayView<double,2> lonlat( mesh->nodes().lonlat() );
  size_t nnodes = mesh->nodes().size();
  size_t nlev = 137;

  FieldSet fields;
  fields.add( fvm.createField<double>("var",nlev) );
  fields.add( fvm.createField<double>("grad",nlev,make_shape(2)) );

  const double deg2rad = M_PI/180.;
  ArrayView<double,2> var( fields["var"] );
  for( size_t jnode=0; jnode< nnodes ; ++jnode )
  {
    const double y  = lonlat(jnode,LAT) * deg2rad ;

    for(size_t jlev = 0; jlev < nlev; ++jlev)
      var(jnode,jlev) = 100.+50.*std::cos(2.*y);
  }
  nabla->gradient(fields["var"],fields["grad"]);

  // output to gmsh
  {
    io::Gmsh().write(*mesh,grid->shortName()+".msh");
    io::Gmsh().write(fields["var"],grid->shortName()+"_fields.msh");
    io::Gmsh().write(fields["grad"],grid->shortName()+"_fields.msh",std::ios::app);
  }
}


} // namespace test
} // namespace atlas
