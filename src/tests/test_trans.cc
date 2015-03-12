/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#define BOOST_TEST_MODULE atlas_test_trans
#include "ecbuild/boost_test_framework.h"

#include <algorithm>    // std::min_element, std::max_element, std::copy

#include "eckit/mpi/mpi.h"
#include "eckit/config/ResourceMgr.h"
#include "atlas/atlas.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/TransPartitioner.h"
#include "atlas/grids/grids.h"
#include "atlas/meshgen/EqualAreaPartitioner.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/LogFormat.h"
#include "atlas/GridDistribution.h"
#include "atlas/io/Gmsh.h"

#include "transi/trans.h"

using namespace eckit;
using namespace atlas;
//using namespace atlas::trans;
using namespace atlas::grids;

struct Fixture   {
       Fixture() {
         trans_init();
         atlas_init(boost::unit_test::framework::master_test_suite().argc,
                    boost::unit_test::framework::master_test_suite().argv);
       }
      ~Fixture() {
         atlas_finalize();
         trans_finalize();
       }
};


void read_rspecg(trans::Trans& trans, std::vector<double>& rspecg, std::vector<int>& nfrom, int &nfld )
{
  eckit::Log::info() << "read_rspecg ...\n";
  nfld = 2;
  if( trans.myproc() == 1 )
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

  eckit::Log::info() << "read_rspecg ... done" << std::endl;
}

//void write_rgpg(Trans* trans, double* rgpg[], int nfld )
//{
//  int jfld;
//  if( trans->myproc == 1 ) printf("write_rgpg ...\n");
//  for( jfld=0; jfld<nfld; ++jfld )
//  {
//    // output global field rgpg[jfld]
//  }
//  if( trans->myproc == 1 ) printf("write_rgpg ... done\n");
//}

namespace atlas {
namespace trans {
class DistSpec
{
public:
  DistSpec(Trans& trans)
  {
    distspec_ = new_distspec(trans);
  }

  DistSpec& nfld(int _nfld) { distspec_.nfld = _nfld; return *this; }

  DistSpec& nfrom (const int*    _nfrom)  { distspec_.nfrom  = const_cast<int*>(_nfrom);     return *this; }
  DistSpec& rspecg(const double* _rspecg) { distspec_.rspecg = const_cast<double*>(_rspecg); return *this; }
  DistSpec& rspec (      double* _rspec)  { distspec_.rspec = _rspec;                        return *this; }

  DistSpec& nfrom (const std::vector<int>&    _nfrom)  { return nfrom(_nfrom.size()?_nfrom.data():NULL); }
  DistSpec& rspecg(const std::vector<double>& _rspecg) { return rspecg(_rspecg.size()?_rspecg.data():NULL); }
  DistSpec& rspec (      std::vector<double>& _rspec)
  { ASSERT(_rspec.size()); return rspec(_rspec.data()); }

  void operator()() { ::trans_distspec(&distspec_); }

private:
  struct DistSpec_t distspec_;
};
}
}

BOOST_GLOBAL_FIXTURE( Fixture )

BOOST_AUTO_TEST_CASE( test_trans_distribution_matches_atlas )
{
  // Create grid and trans object
  ReducedGrid::Ptr g ( ReducedGrid::create( "rgg.N80" ) );

  BOOST_CHECK_EQUAL( g->nlat() , 160 );

  trans::Trans trans( *g );

  BOOST_CHECK_EQUAL( trans.nsmax() , 159 );

  trans::TransPartitioner partitioner(*g,trans);
  GridDistribution distribution( partitioner );

  // -------------- do checks -------------- //
  BOOST_CHECK_EQUAL( trans.nproc(),  eckit::mpi::size() );
  BOOST_CHECK_EQUAL( trans.myproc(), eckit::mpi::rank()+1 );

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

    for( int j=0; j<g->npts(); ++j )
      ++npts[distribution.partition(j)];

    BOOST_CHECK_EQUAL( trans.ngptotg(), g->npts() );
    BOOST_CHECK_EQUAL( trans.ngptot(),  npts[eckit::mpi::rank()] );
    BOOST_CHECK_EQUAL( trans.ngptotmx(), *std::max_element(npts.begin(),npts.end()) );

    for( int j=0; j<partitioner.nb_bands(); ++j )
      BOOST_CHECK_EQUAL( trans.n_regions()[j] , partitioner.nb_regions(j) );
  }

}

BOOST_AUTO_TEST_CASE( test_trans_partitioner )
{
  // Create grid and trans object
  ReducedGrid::Ptr g ( ReducedGrid::create( "rgg.N80" ) );

  trans::Trans trans( *g, 0 );

  BOOST_CHECK_EQUAL( trans.nsmax() , 0 );
  BOOST_CHECK_EQUAL( trans.ngptotg() , g->npts() );

}

BOOST_AUTO_TEST_CASE( test_distspec )
{
  ReducedGrid::Ptr g ( ReducedGrid::create( "rgg4.N80" ) );
  eckit::ResourceMgr::instance().set("atlas.meshgen.angle","0");
  meshgen::ReducedGridMeshGenerator generate;
  BOOST_CHECKPOINT("mesh generator created");
  trans::Trans trans(*g);
  BOOST_CHECKPOINT("Trans initialized");
  std::vector<double> rspecg;
  std::vector<int   > nfrom;
  int nfld;
  BOOST_CHECKPOINT("Read rspecg");
  read_rspecg(trans,rspecg,nfrom,nfld);

  std::vector<double> rspec(nfld*trans.nspec2());

  BOOST_CHECKPOINT("distspec");

  trans::DistSpec distspec(trans);
  distspec.nfld(nfld).nfrom(nfrom).rspecg(rspecg).rspec(rspec)();
  BOOST_CHECKPOINT("end test_distspec");
}

BOOST_AUTO_TEST_CASE( test_distribution )
{
  ReducedGrid::Ptr g ( ReducedGrid::create( "rgg4.N80" ) );

  BOOST_CHECKPOINT("test_distribution");

  GridDistribution::Ptr d_trans( trans::TransPartitioner(*g).distribution() );
  BOOST_CHECKPOINT("trans distribution created");

  GridDistribution::Ptr d_eqreg( meshgen::EqualAreaPartitioner(*g).distribution() );

  BOOST_CHECKPOINT("eqregions distribution created");

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
  ReducedGrid::Ptr g ( ReducedGrid::create( "rgg4.N80" ) );
  eckit::ResourceMgr::instance().set("atlas.meshgen.angle","0");
  eckit::ResourceMgr::instance().set("atlas.meshgen.triangulate","true");

  meshgen::ReducedGridMeshGenerator generate;
  trans::Trans trans(*g);

  Mesh::Ptr m_default( generate( *g ) );
  Mesh::Ptr m_trans( generate( *g, trans::TransPartitioner(*g).distribution() ) );
  Mesh::Ptr m_eqreg( generate( *g, meshgen::EqualAreaPartitioner(*g).distribution() ) );

  ArrayView<int,1> p_default( m_default->function_space("nodes").field("partition") );
  ArrayView<int,1> p_trans  ( m_trans  ->function_space("nodes").field("partition") );
  ArrayView<int,1> p_eqreg  ( m_eqreg  ->function_space("nodes").field("partition") );

  BOOST_CHECK_EQUAL_COLLECTIONS( p_default.begin(), p_default.end(),
                                 p_trans  .begin(), p_trans  .end() );

  BOOST_CHECK_EQUAL_COLLECTIONS( p_default.begin(), p_default.end(),
                                 p_eqreg  .begin(), p_eqreg  .end() );

  //Mesh::Ptr mesh ( generate(*g, meshgen::EqualAreaPartitioner(*g).distribution() ) );

  io::Gmsh().write(*m_trans,"N16_trans.msh");
}


