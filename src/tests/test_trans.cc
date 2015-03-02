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

#include <algorithm>    // std::min_element, std::max_element

#include "eckit/mpi/mpi.h"
#include "trans_api/trans_api.h"
#include "atlas/atlas.h"
#include "atlas/grids/grids.h"
#include "atlas/meshgen/EqualAreaPartitioner.h"
#include "atlas/LogFormat.h"

using namespace eckit;
using namespace atlas;
using namespace atlas::grids;

struct Fixture   {
       Fixture() { atlas_init();     trans_init();     }
      ~Fixture() { atlas_finalize(); trans_finalize(); }
};

namespace atlas {
namespace trans {

class EqualRegionPartitioner
{
public:

  EqualRegionPartitioner( Trans& trans ) :
    trans_(trans),
    partitioner_(trans.nproc)
  {
  }

  void partition(const Grid& grid, int part[]) const
  {

    trans_inquire(&trans_, "nfrstlat,nlstlat,nptrfrstlat,nsta,nonl" );

    int nlonmax = dynamic_cast<const ReducedGrid*>(&grid)->nlonmax();

    int i(0);
    std::vector<int> iglobal(trans_.ndgl*nlonmax);
    for( int jgl=1; jgl<=trans_.ndgl; ++jgl )
    {
      for( int jl=1; jl<=trans_.nloen[jgl-1]; ++jl )
      {
        ++i;
        iglobal[(jgl-1)*nlonmax+(jl-1)] = i;
      }
    }

    const int stride=trans_.ndgl+trans_.n_regions_NS-1;

    int iproc(0);
    for( int ja=1; ja<=trans_.n_regions_NS; ++ja )
    {
      for( int jb=1; jb<=trans_.n_regions[ja-1]; ++jb )
      {
        for( int jgl=trans_.nfrstlat[ja-1]; jgl<=trans_.nlstlat[ja-1]; ++jgl )
        {
          int igl = trans_.nptrfrstlat[ja-1] + jgl - trans_.nfrstlat[ja-1];

          int idx = (jb-1)*stride+(igl-1);
          for( int jl=trans_.nsta[idx]; jl<=trans_.nsta[idx]+trans_.nonl[idx]-1; ++jl )
          {
            int ind = iglobal[(jgl-1)*nlonmax+(jl-1)];
            part[ind-1] = iproc;
          }
        }
        ++iproc;
      }
    }

/*
    IGLOBAL(:,:)=0
    DO JGL=1,NDGLG
      DO JL=1,NLOENG(JGL)
        I=I+1
        IGLOBAL(JL,JGL) = I
      ENDDO
    ENDDO

    DO JA=1,N_REGIONS_NS
      DO JB=1,N_REGIONS(JA)
        IPROC=NGPSET2PE(JA,JB)
        I=0
        DO JGL=NFRSTLAT(JA),NLSTLAT(JA)
          IGL = NPTRFRSTLAT(JA)+JGL-NFRSTLAT(JA)
          DO JL=NSTA(IGL,JB),NSTA(IGL,JB)+NONL(IGL,JB)-1
            IND=IGLOBAL(JL,JGL)
            YRMP%NGLOBALPROC(IND)=IPROC
            I=I+1
            YRMP%NLOCALINDEX(IND)=I
          ENDDO
        ENDDO
      ENDDO
    ENDDO
*/
    //return partitioner_.partition(grid,part);
  }

  int nb_bands() const
  {
    return partitioner_.nb_bands();
  }

  int nb_regions(int b) const
  {
    return partitioner_.nb_regions(b);
  }

private:
  Trans& trans_;
  meshgen::EqualAreaPartitioner partitioner_;

};

}
}

BOOST_GLOBAL_FIXTURE( Fixture )

BOOST_AUTO_TEST_CASE( test_trans_distribution_matches_atlas )
{
  // Create grid and trans object
  ReducedGrid* g = ReducedGrid::create( "rgg.N80" );

  BOOST_CHECK_EQUAL( g->nlat() , 160 );

  Trans trans;

  trans.ndgl  = g->nlat();
  trans.nloen = const_cast<int*>(g->npts_per_lat().data());

  // Assume Linear Grid
  trans.nsmax=(2*trans.ndgl-1)/2;
  BOOST_CHECK_EQUAL( trans.nsmax , 159 );

  // Register resolution in trans library
  trans_setup(&trans);


  // -------------- do checks -------------- //
  BOOST_CHECK_EQUAL( trans.nproc,  eckit::mpi::size() );
  BOOST_CHECK_EQUAL( trans.myproc, eckit::mpi::rank()+1 );

  if( eckit::mpi::rank() == 0 ) // all tasks do the same, so only one needs to check
  {
    trans::EqualRegionPartitioner partitioner(trans);
    int max_nb_regions_EW(0);
    for( int j=0; j<partitioner.nb_bands(); ++j )
      max_nb_regions_EW = std::max(max_nb_regions_EW, partitioner.nb_regions(j));

    BOOST_CHECK_EQUAL( trans.n_regions_NS, partitioner.nb_bands() );
    BOOST_CHECK_EQUAL( trans.n_regions_EW, max_nb_regions_EW );

    std::vector<int> part(g->npts());
    std::vector<int> npts(eckit::mpi::size(),0);

    partitioner.partition(*g,part.data());
    for( int j=0; j<g->npts(); ++j )
      ++npts[part[j]];

    BOOST_CHECK_EQUAL( trans.ngptotg, g->npts() );
    BOOST_CHECK_EQUAL( trans.ngptot, npts[eckit::mpi::rank()] );
    BOOST_CHECK_EQUAL( trans.ngptotmx, *std::max_element(npts.begin(),npts.end()) );

    //for( int j=0; j<partitioner.nb_bands(); ++j )
    //BOOST_CHECK_EQUAL( trans.n_regions[j] , partitioner.nb_regions(j) );
  }
  trans_delete(&trans);
}




