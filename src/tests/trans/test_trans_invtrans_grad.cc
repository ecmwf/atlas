/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#define BOOST_TEST_MODULE atlas_test_trans_invtrans_grad
#include "ecbuild/boost_test_framework.h"

#include <algorithm>

#include "atlas/atlas.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/GridDistribution.h"
#include "atlas/grid/grids.h"
#include "atlas/grid/partitioners/EqualRegionsPartitioner.h"
#include "atlas/grid/partitioners/TransPartitioner.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/trans/Trans.h"
#include "atlas/util/Constants.h"
#include "atlas/internals/Parameters.h"
#include "transi/trans.h"

#include "tests/AtlasFixture.h"

using namespace eckit;
using namespace atlas::internals;
namespace atlas {
namespace test {

struct AtlasTransFixture : public AtlasFixture {
       AtlasTransFixture() {
         trans_init();
       }

      ~AtlasTransFixture() {
         trans_finalize();
       }
};

/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
static void rotated_flow_magnitude(grid::Structured& grid, double var[], const double& beta)
{
  const double radius = util::Earth::radiusInMeters();
  const double USCAL = 20.;
  const double pvel = USCAL/radius;
  const double deg2rad = M_PI/180.;


  size_t n(0);
  for( size_t jlat=0; jlat<grid.nlat(); ++jlat ) {
    for( size_t jlon=0; jlon<grid.nlon(jlat); ++jlon ) {
      const double x = grid.lon(jlat,jlon) * deg2rad;
      const double y = grid.lat(jlat)      * deg2rad;
      const double Ux =  pvel*(std::cos(beta)+std::tan(y)*std::cos(x)*std::sin(beta))*radius*std::cos(y);
      const double Uy = -pvel*std::sin(x)*std::sin(beta)*radius;
      var[n] = std::sqrt(Ux*Ux+Uy*Uy);
      ++n;
    }
  }
}

/// @brief Compute magnitude of flow with rotation-angle beta
/// (beta=0 --> zonal, beta=pi/2 --> meridional)
void rotated_flow_magnitude(const functionspace::NodeColumns& fs, field::Field& field, const double& beta)
{
  const double radius = util::Earth::radiusInMeters();
  const double USCAL = 20.;
  const double pvel = USCAL/radius;
  const double deg2rad = M_PI/180.;

  array::ArrayView<double,2> lonlat_deg = array::make_view<double,2>(fs.nodes().lonlat());
  array::ArrayView<double,1> var        = array::make_view<double,1>(field);

  size_t nnodes = fs.nodes().size();
  for( size_t jnode=0; jnode<nnodes; ++jnode )
  {
     double x = lonlat_deg(jnode,LON) * deg2rad;
     double y = lonlat_deg(jnode,LAT) * deg2rad;
     double Ux =  pvel*(std::cos(beta)+std::tan(y)*std::cos(x)*std::sin(beta))*radius*std::cos(y);
     double Uy = -pvel*std::sin(x)*std::sin(beta)*radius;
     var(jnode) = std::sqrt(Ux*Ux+Uy*Uy);
  }
}


BOOST_GLOBAL_FIXTURE( AtlasTransFixture );


BOOST_AUTO_TEST_CASE( test_invtrans_ifsStyle )
{
  std::string grid_uid("O80");
  SharedPtr<grid::Structured> g ( grid::Structured::create( grid_uid ) );

  trans::Trans trans(*g, g->N()*2-1);
  BOOST_TEST_CHECKPOINT("Trans initialized");
  std::vector<double> rspecg;
  int nfld = 1;

  std::vector<double> init_gpg(trans.ngptotg());
  std::vector<double> init_gp (trans.ngptot ());
  std::vector<double> init_sp (trans.nspec2 ());
  std::vector<int>    nfrom(nfld,1);
  if( parallel::mpi::comm().rank()==0) {
    double beta = M_PI*0.5;
    rotated_flow_magnitude(*g,init_gpg.data(),beta);
  }
  trans.distgrid(nfld,nfrom.data(),init_gpg.data(),init_gp.data());
  trans.dirtrans(nfld,init_gp.data(),init_sp.data());

  std::vector<double> rgp(3*nfld*trans.ngptot());
  double *no_vorticity(NULL), *no_divergence(NULL);
  int nb_vordiv(0);
  int nb_scalar(nfld);
  trans::TransParameters p;
  p.set_scalar_derivatives(true);
  trans.invtrans( nb_scalar, init_sp.data(), nb_vordiv, no_vorticity, no_divergence, rgp.data(), p );

  std::vector<int>    nto(nfld,1);
  std::vector<double> rgpg(3*nfld*trans.ngptotg());

  trans.gathgrid( nfld, nto.data(),   rgp.data(),    rgpg.data() );

  // Output
  {
    mesh::Mesh::Ptr mesh( mesh::generators::Structured().generate(*g) );
    functionspace::StructuredColumns gp(*g);
    output::Gmsh gmsh(grid_uid+"-grid.msh");
    field::Field::Ptr scalar(
          field::Field::wrap<double>("scalar",rgp.data(),array::make_shape(gp.npts())) );
    field::Field::Ptr scalar_dNS(
          field::Field::wrap<double>("scalar_dNS",rgp.data()+nfld*gp.npts(),array::make_shape(gp.npts())));
    field::Field::Ptr scalar_dEW(
          field::Field::wrap<double>("scalar_dEW",rgp.data()+2*nfld*gp.npts(),array::make_shape(gp.npts())));
    gmsh.write(*mesh);
    gmsh.write(*scalar,gp);
    gmsh.write(*scalar_dEW,gp);
    gmsh.write(*scalar_dNS,gp);
  }
}


BOOST_AUTO_TEST_CASE( test_invtrans_grad )
{
  std::string grid_uid("O48");
  SharedPtr<grid::Structured> g ( grid::Structured::create( grid_uid ) );
  mesh::Mesh::Ptr mesh( mesh::generators::Structured().generate(*g) );
  trans::Trans trans(*g, g->N()*2-1);
  functionspace::NodeColumns gp(*mesh);
  functionspace::Spectral sp(trans);

  field::Field::Ptr scalar   ( gp.createField<double>("scalar") );
  field::Field::Ptr scalar_sp( sp.createField<double>("scalar_sp") );
  field::Field::Ptr grad     ( gp.createField<double>("grad",array::make_shape(2)) );

  // Initial condition
  double beta = M_PI*0.5;
  rotated_flow_magnitude(gp,*scalar,beta);

  // Transform to spectral
  trans.dirtrans(gp,*scalar,sp,*scalar_sp);

  // Inverse transform for gradient
  trans.invtrans_grad(sp,*scalar_sp,gp,*grad);

  gp.haloExchange(*grad);

  // Output
  {
    mesh::Mesh::Ptr mesh( mesh::generators::Structured().generate(*g) );
    functionspace::StructuredColumns gp(*g);
    output::Gmsh gmsh(grid_uid+"-nodes.msh");
    gmsh.write(*mesh);
    gmsh.write(*scalar,gp);
    gmsh.write(*grad,gp);
  }
}

} // namespace test
} // namespace atlas
