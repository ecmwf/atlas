/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <algorithm>
#include <cmath>

#define BOOST_TEST_MODULE TestGather
#include "ecbuild/boost_test_framework.h"
#include "eckit/utils/Translator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/internals/atlas_config.h"
#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/internals/Debug.h"

/// POD: Type to test
typedef double POD;

namespace atlas {
namespace test {

template<typename T, size_t N>
std::vector<T> vec( const T (&list)[N] )
{
  return std::vector<T>(list,list+N);
}

struct MPIFixture {
    MPIFixture()  { eckit::mpi::init(); }
    ~MPIFixture() { eckit::mpi::finalize(); }
};

struct Fixture {
  Fixture()
  {
    int nnodes_c[] = {6, 6, 7}; nb_nodes = vec(nnodes_c);
    Nl = nb_nodes[eckit::mpi::rank()];
    switch( eckit::mpi::rank() )
    {
      case 0:
      {               //./----> extra ghost point with nonstandard gidx
        int part_c[] = {2,0,0,0,1,2 };  part = vec(part_c);
        int ridx_c[] = {4,1,2,3,1,3 };  ridx = vec(ridx_c);
        gidx_t gidx_c[] = {9,1,2,3,4,20};  gidx = vec(gidx_c);
      break;
      }
      case 1:
      {
        int part_c[] = {0,1,1,1,2,2};  part = vec(part_c);
        int ridx_c[] = {3,1,2,3,2,3};  ridx = vec(ridx_c);
        gidx_t gidx_c[] = {3,4,5,6,7,8};  gidx = vec(gidx_c);
        break;
      }
      case 2:
      {
        int part_c[] = {1,1,2,2,2,0,0};  part = vec(part_c);
        int ridx_c[] = {2,3,2,3,4,1,2};  ridx = vec(ridx_c);
        gidx_t gidx_c[] = {5,6,7,8,9,1,2};  gidx = vec(gidx_c);
        break;
      }
    }
    gather_scatter.setup(part.data(),ridx.data(),0,gidx.data(),Nl);
    Ng = eckit::mpi::rank() == 0 ? gather_scatter.glb_dof() : 0;
  }
  parallel::GatherScatter gather_scatter;
  std::vector<int> nb_nodes;
  std::vector<int> part;
  std::vector<int> ridx;
  std::vector<gidx_t> gidx;

  int Nl, Ng;
};


BOOST_GLOBAL_FIXTURE( MPIFixture );

BOOST_FIXTURE_TEST_CASE( test_gather_rank0, Fixture )
{
  std::vector<POD> loc(Nl);
  std::vector<POD> glb(Ng);

  for( int j=0; j<Nl; ++j ) {
    loc[j] = (size_t(part[j]) != eckit::mpi::rank() ? 0 : gidx[j]*10 );
  }

  size_t strides[] = {1};
  size_t extents[] = {1};
  gather_scatter.gather(loc.data(),strides,extents,1,glb.data(),strides,extents,1);

  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { 10, 20, 30, 40, 50, 60, 70, 80, 90 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb.begin(),glb.end(), glb_c,glb_c+Ng);
  }
}

BOOST_FIXTURE_TEST_CASE( test_gather_rank1_deprecated, Fixture )
{
  array::ArrayT<POD> loc(Nl,2);
  array::ArrayT<POD> glb(Ng,2);
  array::ArrayT<POD> glb1(Ng,1);
  array::ArrayT<POD> glb2(Ng,1);
  array::ArrayView<POD,2> locv(loc);
  for( int j=0; j<Nl; ++j ) {
    locv(j,0) = (size_t(part[j]) != eckit::mpi::rank() ? 0 : gidx[j]*10 );
    locv(j,1) = (size_t(part[j]) != eckit::mpi::rank() ? 0 : gidx[j]*100);
  }

  // Gather complete field
  {
  size_t loc_strides[] = {1};
  size_t loc_extents[] = {2};
  size_t glb_strides[] = {1};
  size_t glb_extents[] = {2};
  gather_scatter.gather( loc.data(), loc_strides, loc_extents, 1,
                 glb.data(), glb_strides, glb_extents, 1 );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { 10,100, 20,200, 30,300, 40,400, 50,500, 60,600, 70,700, 80,800, 90,900 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb.data(),glb.data()+2*Ng, glb_c,glb_c+2*Ng);
  }

  // Gather only first component
  {
    size_t loc_strides[] = {2};
    size_t loc_extents[] = {1};
    size_t glb_strides[] = {1};
    size_t glb_extents[] = {1};
    gather_scatter.gather( loc.data(),  loc_strides, loc_extents, 1,
                   glb1.data(), glb_strides, glb_extents, 1 );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb1_c[] = { 10, 20, 30, 40, 50, 60, 70, 80, 90 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb1.data(),glb1.data()+Ng, glb1_c,glb1_c+Ng);
  }

  // Gather only second component
  {
    size_t loc_strides[] = {2};
    size_t loc_extents[] = {1};
    size_t glb_strides[] = {1};
    size_t glb_extents[] = {1};
    gather_scatter.gather( loc.data()+1, loc_strides, loc_extents, 1,
                   glb2.data(),  glb_strides, glb_extents, 1 );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb2_c[] = { 100, 200, 300, 400, 500, 600, 700, 800, 900 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb2.data(),glb2.data()+Ng, glb2_c,glb2_c+Ng);
  }
}

BOOST_FIXTURE_TEST_CASE( test_gather_rank1, Fixture )
{
  array::ArrayT<POD> loc(Nl,2);
  array::ArrayT<POD> glb(Ng,2);
  array::ArrayT<POD> glb1(Ng,1);
  array::ArrayT<POD> glb2(Ng,1);
  array::ArrayView<POD,2> locv(loc);
  for( int j=0; j<Nl; ++j ) {
    locv(j,0) = (size_t(part[j]) != eckit::mpi::rank() ? 0 : gidx[j]*10 );
    locv(j,1) = (size_t(part[j]) != eckit::mpi::rank() ? 0 : gidx[j]*100);
  }

  // Gather complete field
  {
  size_t loc_strides[] = {2,1};
  size_t loc_extents[] = {size_t(Nl), 2};
  size_t loc_rank = 2;
  size_t loc_mpl_idxpos[] = {0};
  size_t loc_mpl_rank = 1;
  size_t glb_strides[] = {2,1};
  size_t glb_extents[] = {size_t(Ng), 2};
  size_t glb_rank = 2;
  size_t glb_mpl_idxpos[] = {0};
  size_t glb_mpl_rank = 1;
  size_t root = 0;
  internals::MPL_ArrayView<POD> lview(loc.data(),loc_strides,loc_extents,loc_rank,loc_mpl_idxpos,loc_mpl_rank);
  internals::MPL_ArrayView<POD> gview(glb.data(),glb_strides,glb_extents,glb_rank,glb_mpl_idxpos,glb_mpl_rank);

  BOOST_CHECK_EQUAL(lview.var_rank(),1);
  BOOST_CHECK_EQUAL(lview.var_stride(0),1);
  BOOST_CHECK_EQUAL(lview.var_shape(0),2);
  BOOST_CHECK_EQUAL(gview.var_rank(),1);
  BOOST_CHECK_EQUAL(gview.var_stride(0),1);
  BOOST_CHECK_EQUAL(gview.var_shape(0),2);

  BOOST_CHECK_EQUAL(lview.mpl_rank(),1);
  BOOST_CHECK_EQUAL(lview.mpl_stride(0),2);
  BOOST_CHECK_EQUAL(lview.mpl_shape(0),Nl);
  BOOST_CHECK_EQUAL(gview.mpl_rank(),1);
  BOOST_CHECK_EQUAL(gview.mpl_stride(0),2);
  BOOST_CHECK_EQUAL(gview.mpl_shape(0),Ng);

  gather_scatter.gather( loc.data(), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                         glb.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                         root );

  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { 10,100, 20,200, 30,300, 40,400, 50,500, 60,600, 70,700, 80,800, 90,900 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb.data(),glb.data()+2*Ng, glb_c,glb_c+2*Ng);
  }

  // Gather only first component
  {
    size_t loc_strides[] = {2,2};
    size_t loc_extents[] = {size_t(Nl), 1};
    size_t loc_rank = 2;
    size_t loc_mpl_idxpos[] = {0};
    size_t loc_mpl_rank = 1;
    size_t glb_strides[] = {1};
    size_t glb_extents[] = {size_t(Ng)};
    size_t glb_rank = 1;
    size_t glb_mpl_idxpos[] = {0};
    size_t glb_mpl_rank = 1;
    size_t root = 0;
    internals::MPL_ArrayView<POD> lview(loc.data(),loc_strides,loc_extents,loc_rank,loc_mpl_idxpos,loc_mpl_rank);
    BOOST_CHECK_EQUAL(lview.var_rank(),1);
    BOOST_CHECK_EQUAL(lview.var_stride(0),2);
    BOOST_CHECK_EQUAL(lview.var_shape(0),1);
    BOOST_CHECK_EQUAL(lview.mpl_rank(),1);
    BOOST_CHECK_EQUAL(lview.mpl_stride(0),2);
    BOOST_CHECK_EQUAL(lview.mpl_shape(0),Nl);

    internals::MPL_ArrayView<POD> gview(glb1.data(),glb_strides,glb_extents,glb_rank,glb_mpl_idxpos,glb_mpl_rank);
    BOOST_CHECK_EQUAL(gview.var_rank(),1);
    BOOST_CHECK_EQUAL(gview.var_stride(0),1);
    BOOST_CHECK_EQUAL(gview.var_shape(0),1);
    BOOST_CHECK_EQUAL(gview.mpl_rank(),1);
    BOOST_CHECK_EQUAL(gview.mpl_stride(0),1);
    BOOST_CHECK_EQUAL(gview.mpl_shape(0),Ng);

    gather_scatter.gather( loc.data(),  loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                           glb1.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                           root );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb1_c[] = { 10, 20, 30, 40, 50, 60, 70, 80, 90 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb1.data(),glb1.data()+Ng, glb1_c,glb1_c+Ng);
  }

  // Gather only second component
  {
    size_t loc_strides[] = {2,2};
    size_t loc_extents[] = {size_t(Nl), 1};
    size_t loc_rank = 2;
    size_t loc_mpl_idxpos[] = {0};
    size_t loc_mpl_rank = 1;
    size_t glb_strides[] = {1};
    size_t glb_extents[] = {size_t(Ng)};
    size_t glb_rank = 1;
    size_t glb_mpl_idxpos[] = {0};
    size_t glb_mpl_rank = 1;
    size_t root = 0;
    gather_scatter.gather( loc.data()+1,  loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                           glb2.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                           root );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb2_c[] = { 100, 200, 300, 400, 500, 600, 700, 800, 900 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb2.data(),glb2.data()+Ng, glb2_c,glb2_c+Ng);
  }
}



BOOST_FIXTURE_TEST_CASE( test_gather_rank2, Fixture )
{
  array::ArrayT<POD> loc(Nl,3,2);
  array::ArrayT<POD> glb(Ng,3,2);
  array::ArrayT<POD> glbx1(Ng,3);
  array::ArrayT<POD> glbx2(Ng,3);
  array::ArrayT<POD> glb1x(Ng,2);
  array::ArrayT<POD> glb2x(Ng,2);
  array::ArrayT<POD> glb32(Ng);

  array::ArrayView<POD,3> locv(loc);
  for(int p = 0; p < Nl; ++p)
  {
    for(int i = 0; i < 3; ++i)
    {
      locv(p,i,0) = (size_t(part[p]) != eckit::mpi::rank() ? 0 : -gidx[p]*std::pow(10,i) );
      locv(p,i,1) = (size_t(part[p]) != eckit::mpi::rank() ? 0 :  gidx[p]*std::pow(10,i) );
    }
  }

  // Gather complete field
  {
    size_t loc_strides[] = {6,2,1};
    size_t loc_extents[] = {size_t(Nl), 3, 2};
    size_t loc_rank = 3;
    size_t loc_mpl_idxpos[] = {0};
    size_t loc_mpl_rank = 1;
    size_t glb_strides[] = {6,2,1};
    size_t glb_extents[] = {size_t(Ng), 3, 2};
    size_t glb_rank = 3;
    size_t glb_mpl_idxpos[] = {0};
    size_t glb_mpl_rank = 1;
    size_t root = 0;
    gather_scatter.gather( loc.data(), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                           glb.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                           root );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { -1,1, -10,10, -100,100,
                    -2,2, -20,20, -200,200,
                    -3,3, -30,30, -300,300,
                    -4,4, -40,40, -400,400,
                    -5,5, -50,50, -500,500,
                    -6,6, -60,60, -600,600,
                    -7,7, -70,70, -700,700,
                    -8,8, -80,80, -800,800,
                    -9,9, -90,90, -900,900 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb.data(),glb.data()+6*Ng, glb_c,glb_c+6*Ng);
  }

  // Gather var 1
  {
    size_t loc_strides[] = {6,2,2};
    size_t loc_extents[] = {size_t(Nl), 3, 1};
    size_t loc_rank = 3;
    size_t loc_mpl_idxpos[] = {0};
    size_t loc_mpl_rank = 1;
    size_t glb_strides[] = {6,1};
    size_t glb_extents[] = {size_t(Ng), 3};
    size_t glb_rank = 2;
    size_t glb_mpl_idxpos[] = {0};
    size_t glb_mpl_rank = 1;
    size_t root = 0;
    gather_scatter.gather( &locv(0,0,0), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                           glbx1.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                           root );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { -1, -10, -100,
                    -2, -20, -200,
                    -3, -30, -300,
                    -4, -40, -400,
                    -5, -50, -500,
                    -6, -60, -600,
                    -7, -70, -700,
                    -8, -80, -800,
                    -9, -90, -900 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glbx1.data(),glbx1.data()+3*Ng, glb_c,glb_c+3*Ng);
  }

  // Gather var 2
  {
    size_t loc_strides[] = {6,2,2};
    size_t loc_extents[] = {size_t(Nl), 3, 1};
    size_t loc_rank = 3;
    size_t loc_mpl_idxpos[] = {0};
    size_t loc_mpl_rank = 1;
    size_t glb_strides[] = {6,1};
    size_t glb_extents[] = {size_t(Ng), 3};
    size_t glb_rank = 2;
    size_t glb_mpl_idxpos[] = {0};
    size_t glb_mpl_rank = 1;
    size_t root = 0;
    gather_scatter.gather( &locv(0,0,1), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                           glbx2.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                           root );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { 1, 10, 100,
                    2, 20, 200,
                    3, 30, 300,
                    4, 40, 400,
                    5, 50, 500,
                    6, 60, 600,
                    7, 70, 700,
                    8, 80, 800,
                    9, 90, 900 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glbx2.data(),glbx2.data()+3*Ng, glb_c,glb_c+3*Ng);
  }

  // Gather lev 1
  {
    size_t loc_strides[] = {6,6,1};
    size_t loc_extents[] = {size_t(Nl), 1, 2};
    size_t loc_rank = 3;
    size_t loc_mpl_idxpos[] = {0};
    size_t loc_mpl_rank = 1;
    size_t glb_strides[] = {2,1};
    size_t glb_extents[] = {size_t(Ng), 2};
    size_t glb_rank = 2;
    size_t glb_mpl_idxpos[] = {0};
    size_t glb_mpl_rank = 1;
    size_t root = 0;

    gather_scatter.gather( &locv(0,0,0), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                           glb1x.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                           root );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { -1,1,
                    -2,2,
                    -3,3,
                    -4,4,
                    -5,5,
                    -6,6,
                    -7,7,
                    -8,8,
                    -9,9, };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb1x.data(),glb1x.data()+2*Ng, glb_c,glb_c+2*Ng);
  }

  // Gather lev 2
  {
    size_t loc_strides[] = {6,6,1};
    size_t loc_extents[] = {size_t(Nl), 1, 2};
    size_t loc_rank = 3;
    size_t loc_mpl_idxpos[] = {0};
    size_t loc_mpl_rank = 1;
    size_t glb_strides[] = {2,1};
    size_t glb_extents[] = {size_t(Ng), 2};
    size_t glb_rank = 2;
    size_t glb_mpl_idxpos[] = {0};
    size_t glb_mpl_rank = 1;
    size_t root = 0;
    gather_scatter.gather( &locv(0,1,0), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                           glb2x.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                           root );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { -10,10,
                    -20,20,
                    -30,30,
                    -40,40,
                    -50,50,
                    -60,60,
                    -70,70,
                    -80,80,
                    -90,90, };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb2x.data(),glb2x.data()+2*Ng, glb_c,glb_c+2*Ng);
  }

  // Gather lev 3 var 2
  {
    size_t loc_strides[] = {6,6,2};
    size_t loc_extents[] = {size_t(Nl), 1, 1};
    size_t loc_rank = 3;
    size_t loc_mpl_idxpos[] = {0};
    size_t loc_mpl_rank = 1;
    size_t glb_strides[] = {1};
    size_t glb_extents[] = {size_t(Ng)};
    size_t glb_rank = 1;
    size_t glb_mpl_idxpos[] = {0};
    size_t glb_mpl_rank = 1;
    size_t root = 0;
    gather_scatter.gather( &locv(0,2,1), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                           glb32.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                           root );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { 100,
                    200,
                    300,
                    400,
                    500,
                    600,
                    700,
                    800,
                    900 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb32.data(),glb32.data()+Ng, glb_c,glb_c+Ng);
  }

}


BOOST_FIXTURE_TEST_CASE( test_gather_rank0_ArrayView, Fixture )
{
  array::ArrayT<POD> loc(Nl);
  array::ArrayT<POD> glb(Ng);

  array::ArrayView<POD,1> locv(loc);
  array::ArrayView<POD,1> glbv(glb);
  for(int p = 0; p < Nl; ++p)
  {
    locv(p) = (size_t(part[p]) != eckit::mpi::rank() ? 0 :  gidx[p]*10 );
  }

  // Gather complete field
  {
    gather_scatter.gather( locv, glbv );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { 10,
                    20,
                    30,
                    40,
                    50,
                    60,
                    70,
                    80,
                    90 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb.data(),glb.data()+Ng, glb_c,glb_c+Ng);
  }

}

BOOST_FIXTURE_TEST_CASE( test_gather_rank1_ArrayView, Fixture )
{
  array::ArrayT<POD> loc(Nl,2);
  array::ArrayT<POD> glb(Ng,2);

  array::ArrayView<POD,2> locv(loc);
  array::ArrayView<POD,2> glbv(glb);
  for(int p = 0; p < Nl; ++p)
  {
    locv(p,0) = (size_t(part[p]) != eckit::mpi::rank() ? 0 : -gidx[p]*10 );
    locv(p,1) = (size_t(part[p]) != eckit::mpi::rank() ? 0 :  gidx[p]*10 );
  }

  // Gather complete field
  {
    gather_scatter.gather( locv, glbv );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { -10,10,
                    -20,20,
                    -30,30,
                    -40,40,
                    -50,50,
                    -60,60,
                    -70,70,
                    -80,80,
                    -90,90 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb.data(),glb.data()+2*Ng, glb_c,glb_c+2*Ng);
  }

}


BOOST_FIXTURE_TEST_CASE( test_gather_rank2_ArrayView, Fixture )
{
  array::ArrayT<POD> loc(Nl,3,2);
  array::ArrayT<POD> glb(Ng,3,2);

  array::ArrayView<POD,3> locv(loc);
  array::ArrayView<POD,3> glbv(glb);
  for(int p = 0; p < Nl; ++p)
  {
    for(int i = 0; i < 3; ++i)
    {
      locv(p,i,0) = (size_t(part[p]) != eckit::mpi::rank() ? 0 : -gidx[p]*std::pow(10,i) );
      locv(p,i,1) = (size_t(part[p]) != eckit::mpi::rank() ? 0 :  gidx[p]*std::pow(10,i) );
    }
  }

  // Gather complete field
  {
    gather_scatter.gather( locv, glbv );
  }
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { -1,1, -10,10, -100,100,
                    -2,2, -20,20, -200,200,
                    -3,3, -30,30, -300,300,
                    -4,4, -40,40, -400,400,
                    -5,5, -50,50, -500,500,
                    -6,6, -60,60, -600,600,
                    -7,7, -70,70, -700,700,
                    -8,8, -80,80, -800,800,
                    -9,9, -90,90, -900,900 };
    BOOST_CHECK_EQUAL_COLLECTIONS(glb.data(),glb.data()+6*Ng, glb_c,glb_c+6*Ng);
  }

}

BOOST_FIXTURE_TEST_CASE( test_scatter_rank2_ArrayView, Fixture )
{
  array::ArrayT<POD> loc(Nl,3,2);
  array::ArrayT<POD> glb(Ng,3,2);

  array::ArrayView<POD,3> locv(loc);
  array::ArrayView<POD,3> glbv(glb);
  if( eckit::mpi::rank() == 0 )
  {
    POD glb_c[] = { -1,1, -10,10, -100,100,
                    -2,2, -20,20, -200,200,
                    -3,3, -30,30, -300,300,
                    -4,4, -40,40, -400,400,
                    -5,5, -50,50, -500,500,
                    -6,6, -60,60, -600,600,
                    -7,7, -70,70, -700,700,
                    -8,8, -80,80, -800,800,
                    -9,9, -90,90, -900,900 };
    glb.assign(glb_c,glb_c+Ng*6);
  }

  POD nan = -1000.;
  locv = nan;

  gather_scatter.scatter( glbv, locv );

  switch( eckit::mpi::rank() )
  {
  case 0: {
    POD loc_c[] = { nan,nan, nan,nan, nan,nan,
                     -1,1,   -10,10, -100,100,
                     -2,2,   -20,20, -200,200,
                     -3,3,   -30,30, -300,300,
                    nan,nan, nan,nan, nan,nan,
                    nan,nan, nan,nan, nan,nan };
    BOOST_CHECK_EQUAL_COLLECTIONS(loc.data(),loc.data()+Nl*6, loc_c,loc_c+Nl*6);
    break; }
  case 1: {
    POD loc_c[] = { nan,nan, nan,nan, nan,nan,
                     -4,4,   -40,40, -400,400,
                     -5,5,   -50,50, -500,500,
                     -6,6,   -60,60, -600,600,
                    nan,nan, nan,nan, nan,nan,
                    nan,nan, nan,nan, nan,nan };
    BOOST_CHECK_EQUAL_COLLECTIONS(loc.data(),loc.data()+Nl*6, loc_c,loc_c+Nl*6);
    break; }
  case 2: {
      POD loc_c[] = { nan,nan, nan,nan, nan,nan,
                      nan,nan, nan,nan, nan,nan,
                       -7,7,   -70,70, -700,700,
                       -8,8,   -80,80, -800,800,
                       -9,9,   -90,90, -900,900,
                      nan,nan, nan,nan, nan,nan,
                      nan,nan, nan,nan, nan,nan };
      BOOST_CHECK_EQUAL_COLLECTIONS(loc.data(),loc.data()+Nl*6, loc_c,loc_c+Nl*6);
      break; }
  }
}

} // namespace test
} // namespace atlas

