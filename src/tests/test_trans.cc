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
#include "atlas/grids/grids.h"
#include "atlas/meshgen/EqualAreaPartitioner.h"
#include "atlas/meshgen/ReducedGridMeshGenerator.h"
#include "atlas/LogFormat.h"
#include "atlas/GridDistribution.h"

#include "trans_api/trans_api.h"


using namespace eckit;
using namespace atlas;
using namespace atlas::grids;

namespace atlas {
namespace trans {

class Trans {
private:
  typedef struct ::Trans_t Trans_t;
public:
  Trans(const ReducedGrid& g)
  {
    int nsmax = (2*g.nlat()-1)/2;
    ctor(g.nlat(),g.npts_per_lat().data(), nsmax);
  }

  Trans(const ReducedGrid& g, const int nsmax )
  {
   ctor(g.nlat(),g.npts_per_lat().data(), nsmax);
  }

  Trans( const std::vector<int>& npts_per_lat, const int nsmax )
  {
    ctor(npts_per_lat.size(),npts_per_lat.data(), nsmax);
  }

  virtual ~Trans()
  {
    ::trans_delete(&trans_);
  }

  operator const Trans_t*() const { return &trans_; }

  int        nproc()        const { return trans_.nproc; }
  int        myproc()       const { return trans_.myproc; }
  int        ndgl()         const { return trans_.ndgl; }
  int        nsmax()        const { return trans_.nsmax; }
  int        ngptot()       const { return trans_.ngptot; }
  int        ngptotg()      const { return trans_.ngptotg; }
  int        ngptotmx()     const { return trans_.ngptotmx; }
  const int* nloen()        const { ASSERT( trans_.nloen     != NULL ); return trans_.nloen; }
  const int* n_regions()    const { ASSERT( trans_.n_regions != NULL ); return trans_.n_regions; }
  int        n_regions_NS() const { return trans_.n_regions_NS; }
  int        n_regions_EW() const { return trans_.n_regions_EW; }
  const int* nfrstlat()     const { if( trans_.nfrstlat    == NULL ) ::trans_inquire(&trans_,"nfrstlat");    return trans_.nfrstlat; }
  const int* nlstlat()      const { if( trans_.nlstlat     == NULL ) ::trans_inquire(&trans_,"nlstlat");     return trans_.nlstlat; }
  const int* nptrfrstlat()  const { if( trans_.nptrfrstlat == NULL ) ::trans_inquire(&trans_,"nptrfrstlat"); return trans_.nptrfrstlat; }
  const int* nsta()         const { if( trans_.nsta        == NULL ) ::trans_inquire(&trans_,"nsta");        return trans_.nsta; }
  const int* nonl()         const { if( trans_.nonl        == NULL ) ::trans_inquire(&trans_,"nonl");        return trans_.nonl; }

private:

  void ctor(const int ndgl, const int nloen[], int nsmax)
  {
    trans_.ndgl  = ndgl;
    trans_.nloen = new int[trans_.ndgl];
    std::copy(nloen,nloen+ndgl,trans_.nloen);
    trans_.nsmax = nsmax;
    trans_setup(&trans_);
  }

private:
  mutable Trans_t trans_;
};

class TransPartitioner: public Partitioner
{
public:

  TransPartitioner( Trans& trans ) :
    Partitioner(),
    t_(trans)
  {
    Partitioner::set_nb_partition( t_.nproc() );
  }

  virtual void partition(const Grid& grid, int part[]) const
  {
    if( dynamic_cast<const ReducedGrid*>(&grid) == NULL )
      throw BadCast("Grid is not a ReducedGrid type. Cannot partition using IFS trans",Here());

    int nlonmax = dynamic_cast<const ReducedGrid*>(&grid)->nlonmax();

    int i(0);
    std::vector<int> iglobal(t_.ndgl()*nlonmax);
    for( int jgl=1; jgl<=t_.ndgl(); ++jgl )
    {
      for( int jl=1; jl<=t_.nloen()[jgl-1]; ++jl )
      {
        ++i;
        iglobal[(jgl-1)*nlonmax+(jl-1)] = i;
      }
    }

    const int stride=t_.ndgl()+t_.n_regions_NS()-1;

    int iproc(0);
    for( int ja=1; ja<=t_.n_regions_NS(); ++ja )
    {
      for( int jb=1; jb<=t_.n_regions()[ja-1]; ++jb )
      {
        for( int jgl=t_.nfrstlat()[ja-1]; jgl<=t_.nlstlat()[ja-1]; ++jgl )
        {
          int igl = t_.nptrfrstlat()[ja-1] + jgl - t_.nfrstlat()[ja-1];

          int idx = (jb-1)*stride+(igl-1);
          for( int jl=t_.nsta()[idx]; jl<=t_.nsta()[idx]+t_.nonl()[idx]-1; ++jl )
          {
            int ind = iglobal[(jgl-1)*nlonmax+(jl-1)];
            part[ind-1] = iproc;
          }
        }
        ++iproc;
      }
    }
  }

  int nb_bands() const
  {
    return t_.n_regions_NS();
  }

  int nb_regions(int b) const
  {
    return t_.n_regions()[b];
  }

private:
  Trans& t_;
};

}
}




using namespace eckit;
using namespace atlas;
//using namespace atlas::trans;
using namespace atlas::grids;

struct Fixture   {
       Fixture() { trans_init(); atlas_init();          }
      ~Fixture() { atlas_finalize(); trans_finalize();  }
};






BOOST_GLOBAL_FIXTURE( Fixture )

BOOST_AUTO_TEST_CASE( test_trans_distribution_matches_atlas )
{
  // Create grid and trans object
  ReducedGrid::Ptr g ( ReducedGrid::create( "rgg.N80" ) );

  BOOST_CHECK_EQUAL( g->nlat() , 160 );

  trans::Trans trans( *g );

  BOOST_CHECK_EQUAL( trans.nsmax() , 159 );

  trans::TransPartitioner partitioner(trans);
  GridDistribution distribution( *g, partitioner );

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


BOOST_AUTO_TEST_CASE( test_generate_mesh )
{
  ReducedGrid::Ptr g ( ReducedGrid::create( "rgg4.N16" ) );
  eckit::ResourceMgr::instance().set("atlas.meshgen.angle","0");
  meshgen::ReducedGridMeshGenerator generate;
  trans::Trans trans(*g);
  //generate(*g,trans::TransPartitioner(trans));
}


