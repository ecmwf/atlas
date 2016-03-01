/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#define BOOST_TEST_MODULE atlas_test_trans
#include "ecbuild/boost_test_framework.h"

#include <algorithm>

#include "eckit/config/ResourceMgr.h"
#include "atlas/util/parallel/mpi/mpi.h"
#include "atlas/atlas.h"
#include "atlas/numerics/trans/Trans.h"
#include "atlas/grid/partitioners/TransPartitioner.h"
#include "atlas/grid/grids.h"
#include "atlas/grid/partitioners/EqualRegionsPartitioner.h"
#include "atlas/mesh/generators/ReducedGridMeshGenerator.h"
#include "atlas/util/runtime/LogFormat.h"
#include "atlas/grid/GridDistribution.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/field/FieldSet.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/ReducedGridColumns.h"

#include "transi/trans.h"

using namespace eckit;
using namespace atlas::grid;

namespace atlas {
namespace test {

struct Fixture   {
       Fixture() {

         atlas_init(boost::unit_test::framework::master_test_suite().argc,
                    boost::unit_test::framework::master_test_suite().argv);
         trans_init();

       }
      ~Fixture() {
         trans_finalize();
         atlas_finalize();
       }
};


void read_rspecg(numerics::trans::Trans& trans, std::vector<double>& rspecg, std::vector<int>& nfrom, int &nfld )
{
  Log::info() << "read_rspecg ...\n";
  nfld = 2;
  if( trans.myproc(0) == 0 )
  {
    rspecg.resize(nfld*trans.nspec2g());
    for( int i=0; i<trans.nspec2g(); ++i )
    {
      rspecg[i*nfld + 0] = (i==0 ? 1. : 0.); // scalar field 1
      rspecg[i*nfld + 1] = (i==0 ? 2. : 0.); // scalar field 2
    }
  }
  nfrom.resize(nfld);
  for (int jfld=0; jfld<nfld; ++jfld)
    nfrom[jfld] = 1;

  Log::info() << "read_rspecg ... done" << std::endl;
}


BOOST_GLOBAL_FIXTURE( Fixture );

BOOST_AUTO_TEST_CASE( test_trans_distribution_matches_atlas )
{
  BOOST_CHECK( grid::partitioners::PartitionerFactory::has("Trans") );


  // Create grid and trans object
  ReducedGrid::Ptr g ( ReducedGrid::create( "N80" ) );

  BOOST_CHECK_EQUAL( g->nlat() , 160 );

  numerics::trans::Trans trans( *g );

  BOOST_CHECK_EQUAL( trans.nsmax() , 0 );

  grid::partitioners::TransPartitioner partitioner(*g,trans);
  GridDistribution distribution( partitioner );

  // -------------- do checks -------------- //
  BOOST_CHECK_EQUAL( trans.nproc(),  eckit::mpi::size() );
  BOOST_CHECK_EQUAL( trans.myproc(0), eckit::mpi::rank() );


  if( eckit::mpi::rank() == 0 ) // all tasks do the same, so only one needs to check
  {
    int max_nb_regions_EW(0);
    for( int j=0; j<partitioner.nb_bands(); ++j )
      max_nb_regions_EW = std::max(max_nb_regions_EW, partitioner.nb_regions(j));

    BOOST_CHECK_EQUAL( trans.n_regions_NS(), partitioner.nb_bands() );
    BOOST_CHECK_EQUAL( trans.n_regions_EW(), max_nb_regions_EW );

    BOOST_CHECK_EQUAL( distribution.nb_partitions(), eckit::mpi::size() );
    BOOST_CHECK_EQUAL( distribution.partition().size(), g->npts() );

    std::vector<int> npts(distribution.nb_partitions(),0);

    for(size_t j = 0; j < g->npts(); ++j)
      ++npts[distribution.partition(j)];

    BOOST_CHECK_EQUAL( trans.ngptotg(), g->npts() );
    BOOST_CHECK_EQUAL( trans.ngptot(),  npts[eckit::mpi::rank()] );
    BOOST_CHECK_EQUAL( trans.ngptotmx(), *std::max_element(npts.begin(),npts.end()) );

    array::ArrayView<int,1> n_regions ( trans.n_regions() ) ;
    for( int j=0; j<partitioner.nb_bands(); ++j )
      BOOST_CHECK_EQUAL( n_regions[j] , partitioner.nb_regions(j) );
  }
}


BOOST_AUTO_TEST_CASE( test_trans_partitioner )
{
  BOOST_TEST_CHECKPOINT("test_trans_partitioner");
  // Create grid and trans object
  ReducedGrid::Ptr g ( ReducedGrid::create( "N80" ) );

  numerics::trans::Trans trans( *g, 0 );

  BOOST_CHECK_EQUAL( trans.nsmax() , 0 );
  BOOST_CHECK_EQUAL( trans.ngptotg() , g->npts() );
}

BOOST_AUTO_TEST_CASE( test_trans_options )
{
  numerics::trans::Trans::Options opts;
  opts.set_fft(numerics::trans::FFTW);
  opts.set_split_latitudes(false);
  opts.set_read("readfile");

  Log::info() << "trans_opts = " << opts << std::endl;
}


BOOST_AUTO_TEST_CASE( test_distspec )
{
  ReducedGrid::Ptr g ( ReducedGrid::create( "O80" ) );
  eckit::ResourceMgr::instance().set("atlas.meshgen.angle","0");
  mesh::generators::ReducedGridMeshGenerator generate;
  BOOST_TEST_CHECKPOINT("mesh generator created");
  //numerics::trans::Trans trans(*g, 159 );

  numerics::trans::Trans::Options p;
  if( eckit::mpi::size() == 1 )
    p.set_write("cached_legendre_coeffs");
  p.set_flt(false);
  numerics::trans::Trans trans(400, 159, p);
  BOOST_TEST_CHECKPOINT("Trans initialized");
  std::vector<double> rspecg;
  std::vector<int   > nfrom;
  int nfld;
  BOOST_TEST_CHECKPOINT("Read rspecg");
  read_rspecg(trans,rspecg,nfrom,nfld);


  std::vector<double> rspec(nfld*trans.nspec2());
  std::vector<int> nto(nfld,1);
  std::vector<double> rgp(nfld*trans.ngptot());
  std::vector<double> rgpg(nfld*trans.ngptotg());

  trans.distspec( nfld, nfrom.data(), rspecg.data(), rspec.data() );
  trans.invtrans( nfld, rspec.data(), rgp.data() );
  trans.gathgrid( nfld, nto.data(),   rgp.data(),    rgpg.data() );

  BOOST_TEST_CHECKPOINT("end test_distspec");
}

BOOST_AUTO_TEST_CASE( test_distribution )
{
  ReducedGrid::Ptr g ( ReducedGrid::create( "O80" ) );

  BOOST_TEST_CHECKPOINT("test_distribution");

  GridDistribution::Ptr d_trans( grid::partitioners::TransPartitioner(*g).distribution() );
  BOOST_TEST_CHECKPOINT("trans distribution created");

  GridDistribution::Ptr d_eqreg( grid::partitioners::EqualRegionsPartitioner(*g).distribution() );

  BOOST_TEST_CHECKPOINT("eqregions distribution created");

  if( eckit::mpi::rank() == 0 )
  {
    BOOST_CHECK_EQUAL( d_trans->nb_partitions(), d_eqreg->nb_partitions() );
    BOOST_CHECK_EQUAL( d_trans->max_pts(), d_eqreg->max_pts() );
    BOOST_CHECK_EQUAL( d_trans->min_pts(), d_eqreg->min_pts() );

    BOOST_CHECK_EQUAL_COLLECTIONS( d_trans->nb_pts().begin(), d_trans->nb_pts().end(),
                                   d_eqreg->nb_pts().begin(), d_eqreg->nb_pts().end() );
  }

}


BOOST_AUTO_TEST_CASE( test_generate_mesh )
{
  BOOST_TEST_CHECKPOINT("test_generate_mesh");
  ReducedGrid::Ptr g ( ReducedGrid::create( "O80" ) );
  eckit::ResourceMgr::instance().set("atlas.meshgen.angle","0");
  eckit::ResourceMgr::instance().set("atlas.meshgen.triangulate","true");

  mesh::generators::ReducedGridMeshGenerator generate;
  numerics::trans::Trans trans(*g);

  mesh::Mesh::Ptr m_default( generate( *g ) );

  BOOST_TEST_CHECKPOINT("trans_distribution");
  GridDistribution::Ptr trans_distribution( grid::partitioners::TransPartitioner(*g).distribution() );
  mesh::Mesh::Ptr m_trans( generate( *g, *trans_distribution ) );

  BOOST_TEST_CHECKPOINT("eqreg_distribution");
  GridDistribution::Ptr eqreg_distribution( grid::partitioners::EqualRegionsPartitioner(*g).distribution() );
  mesh::Mesh::Ptr m_eqreg( generate( *g, *eqreg_distribution ) );

  array::ArrayView<int,1> p_default( m_default->nodes().partition() );
  array::ArrayView<int,1> p_trans  ( m_trans  ->nodes().partition() );
  array::ArrayView<int,1> p_eqreg  ( m_eqreg  ->nodes().partition() );

  BOOST_CHECK_EQUAL_COLLECTIONS( p_default.begin(), p_default.end(),
                                 p_trans  .begin(), p_trans  .end() );

  BOOST_CHECK_EQUAL_COLLECTIONS( p_default.begin(), p_default.end(),
                                 p_eqreg  .begin(), p_eqreg  .end() );

  //mesh::Mesh::Ptr mesh ( generate(*g, mesh::generators::EqualAreaPartitioner(*g).distribution() ) );

  util::io::Gmsh().write(*m_trans,"N16_trans.msh");
}


BOOST_AUTO_TEST_CASE( test_spectral_fields )
{
  BOOST_TEST_CHECKPOINT("test_spectral_fields");

  ReducedGrid::Ptr g ( ReducedGrid::create( "O48" ) );
  eckit::ResourceMgr::instance().set("atlas.meshgen.angle","0");
  eckit::ResourceMgr::instance().set("atlas.meshgen.triangulate","false");

  mesh::generators::ReducedGridMeshGenerator generate;
  mesh::Mesh::Ptr m( generate( *g ) );

  numerics::trans::Trans trans(*g,47);


  SharedPtr<functionspace::NodeColumns> nodal (new functionspace::NodeColumns(*m));
  SharedPtr<functionspace::Spectral> spectral (new functionspace::Spectral(trans));

  SharedPtr<field::Field> spf ( spectral->createField<double>("spf") );
  SharedPtr<field::Field> gpf ( nodal->createField<double>("gpf") );


  BOOST_CHECK_NO_THROW( trans.dirtrans(*nodal,*gpf,*spectral,*spf) );
  BOOST_CHECK_NO_THROW( trans.invtrans(*spectral,*spf,*nodal,*gpf) );

  field::FieldSet gpfields;   gpfields.add(*gpf);
  field::FieldSet spfields;   spfields.add(*spf);

  BOOST_CHECK_NO_THROW( trans.dirtrans(*nodal,gpfields,*spectral,spfields) );
  BOOST_CHECK_NO_THROW( trans.invtrans(*spectral,spfields,*nodal,gpfields) );

  gpfields.add(*gpf);
  BOOST_CHECK_THROW(trans.dirtrans(*nodal,gpfields,*spectral,spfields),eckit::SeriousBug);

}


BOOST_AUTO_TEST_CASE( test_nomesh )
{
  BOOST_TEST_CHECKPOINT("test_spectral_fields");

  SharedPtr<ReducedGrid> g ( ReducedGrid::create( "O48" ) );
  SharedPtr<numerics::trans::Trans> trans ( new numerics::trans::Trans(*g,47) );

  SharedPtr<functionspace::Spectral>    spectral    (new functionspace::Spectral(*trans));
  SharedPtr<functionspace::ReducedGridColumns> gridpoints (new functionspace::ReducedGridColumns(*g));

  SharedPtr<field::Field> spfg ( spectral->createGlobalField<double>("spf") );
  SharedPtr<field::Field> spf  ( spectral->createField<double>("spf") );
  SharedPtr<field::Field> gpf  ( gridpoints->createField<double>("gpf") );
  SharedPtr<field::Field> gpfg ( gridpoints->createGlobalField<double>("gpf") );

  array::ArrayView<double,1> spg (*spfg);
  spg = 0.;
  spg(0) = 4.;

  BOOST_CHECK_NO_THROW( spectral->scatter(*spfg,*spf) );

  if( eckit::mpi::rank() == 0 ) {
    array::ArrayView<double,1> sp (*spf);
    BOOST_CHECK_CLOSE( sp(0), 4., 0.001 );
    for( size_t jp=0; jp<sp.size(); ++jp ) {
      Log::debug(2) << "sp("<< jp << ")   :   " << sp(jp) << std::endl;
    }
  }

  BOOST_CHECK_NO_THROW( trans->invtrans(*spf,*gpf) );

  BOOST_CHECK_NO_THROW( gridpoints->gather(*gpf,*gpfg) );

  if( eckit::mpi::rank() == 0 ) {
    array::ArrayView<double,1> gpg (*gpfg);
    for( size_t jp=0; jp<gpg.size(); ++jp ) {
      BOOST_CHECK_CLOSE( gpg(jp), 4., 0.001 );
      Log::debug(3) << "gpg("<<jp << ")   :   " << gpg(jp) << std::endl;
    }
  }

  BOOST_CHECK_NO_THROW( gridpoints->scatter(*gpfg,*gpf) );

  BOOST_CHECK_NO_THROW( trans->dirtrans(*gpf,*spf) );

  BOOST_CHECK_NO_THROW( spectral->gather(*spf,*spfg) );

  if( eckit::mpi::rank() == 0 ) {
    BOOST_CHECK_CLOSE( spg(0), 4., 0.001 );
  }
}


} // namespace test
} // namespace atlas
