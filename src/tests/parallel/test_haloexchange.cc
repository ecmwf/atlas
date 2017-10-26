/*
 * (C) Copyright 1996-2017 ECMWF.
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

#include "eckit/memory/SharedPtr.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/library/config.h"
#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/parallel/HaloExchange.h"


#include "tests/AtlasTestEnvironment.h"
#include "eckit/testing/Test.h"

using namespace eckit::testing;

/// POD: Type to test
typedef double POD;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

template<typename T, size_t N>
std::vector<T> vec( const T (&list)[N] )
{
  return std::vector<T>(list,list+N);
}

struct Fixture {
  Fixture()
  {
    int nnodes_c[] = {5, 6, 7}; nb_nodes = vec(nnodes_c);
    N = nb_nodes[parallel::mpi::comm().rank()];
    switch( parallel::mpi::comm().rank() )
    {
      case 0:
      {
        int part_c[] = {2,0,0,0,1};  part = vec(part_c);
        int ridx_c[] = {4,1,2,3,1};  ridx = vec(ridx_c);
        POD gidx_c[] = {0,1,2,3,0};  gidx = vec(gidx_c);
        break;
      }
      case 1:
      {
        int part_c[] = {0,1,1,1,2,2};  part = vec(part_c);
        int ridx_c[] = {3,1,2,3,2,3};  ridx = vec(ridx_c);
        POD gidx_c[] = {0,4,5,6,0,0};  gidx = vec(gidx_c);
        break;
      }
      case 2:
      {
        int part_c[] = {1,1,2,2,2,0,0};  part = vec(part_c);
        int ridx_c[] = {2,3,2,3,4,1,2};  ridx = vec(ridx_c);
        POD gidx_c[] = {0,0,7,8,9,0,0};  gidx = vec(gidx_c);
        break;
      }
    }
    halo_exchange.setup(part.data(),ridx.data(),0,N);
  }
  parallel::HaloExchange halo_exchange;
  std::vector<int> nb_nodes;
  std::vector<int> part;
  std::vector<int> ridx;
  std::vector<POD>   gidx;

  int N;
};

//-----------------------------------------------------------------------------

CASE("test_haloexchange") {
  SETUP("Fixture") {
    Fixture f;



    SECTION( "test_rank0" )
    {
//      size_t strides[] = {1};
//      size_t shape[] = {1};
//      f.halo_exchange.execute(f.gidx.data(),strides,shape,1);

//      switch( parallel::mpi::comm().rank() )
//      {
//        case 0: { POD gidx_c[] = { 9, 1, 2, 3, 4};
//          EXPECT(f.gidx == make_view(gidx_c,gidx_c+f.N)); break; }
//        case 1: { POD gidx_c[] = { 3, 4, 5, 6, 7, 8};
//          EXPECT(f.gidx == make_view(gidx_c,gidx_c+f.N)); break; }
//        case 2: { POD gidx_c[] = { 5, 6, 7, 8, 9, 1, 2};
//          EXPECT(f.gidx == make_view(gidx_c,gidx_c+f.N)); break; }
//      }
    }

    SECTION( "test_rank1" )
    {
      array::ArrayT<POD> arr(f.N,2);
      array::ArrayView<POD,2> arrv = array::make_view<POD,2>(arr);
      for( int j=0; j<f.N; ++j ) {
        arrv(j,0) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j]*10 );
        arrv(j,1) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j]*100);
      }

      size_t strides[] = {1};
      size_t shape[] = {2};
      f.halo_exchange.execute(arrv,strides,shape,1);

      switch( parallel::mpi::comm().rank() )
      {
        case 0: { POD arr_c[] = { 90,900, 10,100, 20,200, 30,300, 40,400 };
          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
        case 1: { POD arr_c[] = { 30,300, 40,400, 50,500, 60,600, 70,700, 80,800};
          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
        case 2: { POD arr_c[] = { 50,500, 60,600, 70,700, 80,800, 90,900, 10,100, 20,200};
          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
      }
    }

    SECTION( "test_rank1_strided_v1" )
    {
      array::ArrayT<POD> arr(f.N,2);
      array::ArrayView<POD,2> arrv = array::make_view<POD,2>(arr);
      for( int j=0; j<f.N; ++j ) {
        arrv(j,0) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j]*10 );
        arrv(j,1) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j]*100);
      }

      size_t strides[] = {2};
      size_t shape[] = {1};
      f.halo_exchange.execute(arrv,strides,shape,1);

      switch( parallel::mpi::comm().rank() )
      {
        case 0: { POD arr_c[] = { 90,0, 10,100, 20,200, 30,300, 40,0 };
          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
        case 1: { POD arr_c[] = { 30,0, 40,400, 50,500, 60,600, 70,0, 80,0};
          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
        case 2: { POD arr_c[] = { 50,0, 60,0, 70,700, 80,800, 90,900, 10,0, 20,0};
          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
      }
    }

    SECTION( "test_rank1_strided_v2" )
    {
//      array::ArrayT<POD> arr(f.N,2);
//      array::ArrayView<POD,2> arrv = array::make_view<POD,2>(arr);
//      for( int j=0; j<f.N; ++j ) {
//        arrv(j,0) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j]*10 );
//        arrv(j,1) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j]*100);
//      }

//      size_t strides[] = {2};
//      size_t shape[] = {1};
//      f.halo_exchange.execute(&arrv(0,1),strides,shape,1);

//      switch( parallel::mpi::comm().rank() )
//      {
//        case 0: { POD arr_c[] = { 0,900, 10,100, 20,200, 30,300, 0,400 };
//          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
//        case 1: { POD arr_c[] = { 0,300, 40,400, 50,500, 60,600, 0,700, 0,800};
//          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
//        case 2: { POD arr_c[] = { 0,500, 0,600, 70,700, 80,800, 90,900, 0,100, 0,200};
//          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
//      }
    }

    SECTION( "test_rank2" )
    {
      array::ArrayT<POD> arr(f.N,3,2);
      array::ArrayView<POD,3> arrv = array::make_view<POD,3>(arr);
      for( int p=0; p<f.N; ++p )
      {
        for( size_t i=0; i<3; ++i )
        {
          arrv(p,i,0) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 : -f.gidx[p]*std::pow(10,i) );
          arrv(p,i,1) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 :  f.gidx[p]*std::pow(10,i) );
        }
      }

      size_t strides[] = {2,1};
      size_t shape[] = {3,2};
      f.halo_exchange.execute(arrv,strides,shape,2);

      switch( parallel::mpi::comm().rank() )
      {
        case 0:
        {
          int arr_c[] = { -9,9, -90,90, -900,900,
                          -1,1, -10,10, -100,100,
                          -2,2, -20,20, -200,200,
                          -3,3, -30,30, -300,300,
                          -4,4, -40,40, -400,400};
          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
          break;
        }
        case 1:
        {
          int arr_c[] = { -3,3, -30,30, -300,300,
                          -4,4, -40,40, -400,400,
                          -5,5, -50,50, -500,500,
                          -6,6, -60,60, -600,600,
                          -7,7, -70,70, -700,700,
                          -8,8, -80,80, -800,800};
          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
          break;
        }
        case 2:
        {
          int arr_c[] = { -5,5, -50,50, -500,500,
                          -6,6, -60,60, -600,600,
                          -7,7, -70,70, -700,700,
                          -8,8, -80,80, -800,800,
                          -9,9, -90,90, -900,900,
                          -1,1, -10,10, -100,100,
                          -2,2, -20,20, -200,200};
          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
          break;
        }
      }
    }

    SECTION( "test_rank2_l1" )
    {
      array::ArrayT<POD> arr(f.N,3,2);
      array::ArrayView<POD,3> arrv = array::make_view<POD,3>(arr);
      for( int p=0; p<f.N; ++p )
      {
        for( size_t i=0; i<3; ++i )
        {
          arrv(p,i,0) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 : -f.gidx[p]*std::pow(10,i) );
          arrv(p,i,1) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 :  f.gidx[p]*std::pow(10,i) );
        }
      }

      size_t strides[] = {6,1};
      size_t shape[] = {1,2};
      f.halo_exchange.execute(arrv,strides,shape,2);

      switch( parallel::mpi::comm().rank() )
      {
        case 0:
        {
          POD arr_c[] = { -9,9,   0, 0,    0,  0,  // halo
                          -1,1, -10,10, -100,100,  // core
                          -2,2, -20,20, -200,200,  // core
                          -3,3, -30,30, -300,300,  // core
                          -4,4,   0, 0,    0,  0}; // halo
          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
          break;
        }
        case 1:
        {
          POD arr_c[] = { -3,3,   0, 0,    0,  0,  // halo
                          -4,4, -40,40, -400,400,  // core
                          -5,5, -50,50, -500,500,  // core
                          -6,6, -60,60, -600,600,  // core
                          -7,7,   0, 0,    0,  0,  // halo
                          -8,8,   0, 0,    0,  0}; // halo
          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
          break;
        }
        case 2:
        {
          POD arr_c[] = { -5,5,   0, 0,    0,  0,  // halo
                          -6,6,   0, 0,    0,  0,  // halo
                          -7,7, -70,70, -700,700,  // core
                          -8,8, -80,80, -800,800,  // core
                          -9,9, -90,90, -900,900,  // core
                          -1,1,   0, 0,    0,  0,  // halo
                          -2,2,   0, 0,    0,  0}; // halo
          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
          break;
        }
      }
    }

    SECTION( "test_rank2_l2_v2" )
    {
//      // Test rank 2 halo-exchange
//      array::ArrayT<POD> arr(f.N,3,2);
//      array::ArrayView<POD,3> arrv = array::make_view<POD,3>(arr);
//      for( int p=0; p<f.N; ++p )
//      {
//        for( size_t i=0; i<3; ++i )
//        {
//          arrv(p,i,0) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 : -f.gidx[p]*std::pow(10,i) );
//          arrv(p,i,1) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 :  f.gidx[p]*std::pow(10,i) );
//        }
//      }

//      size_t strides[] = {6,2};
//      size_t shape[] = {1,1};
//      f.halo_exchange.execute(&arrv(0,1,1),strides,shape,2);

//      switch( parallel::mpi::comm().rank() )
//      {
//        case 0:
//        {
//          POD arr_c[] = {  0,0,   0,90,    0,  0,  // halo
//                          -1,1, -10,10, -100,100,  // core
//                          -2,2, -20,20, -200,200,  // core
//                          -3,3, -30,30, -300,300,  // core
//                          0,0,   0,40,    0,  0}; // halo
//          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
//          break;
//        }
//        case 1:
//        {
//          POD arr_c[] = {  0,0,   0,30,    0,  0,  // halo
//                          -4,4, -40,40, -400,400,  // core
//                          -5,5, -50,50, -500,500,  // core
//                          -6,6, -60,60, -600,600,  // core
//                          0,0,   0,70,    0,  0,  // halo
//                          0,0,   0,80,    0,  0}; // halo
//          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
//          break;
//        }
//        case 2:
//        {
//          POD arr_c[] = {  0,0,   0,50,    0,  0,  // halo
//                          0,0,   0,60,    0,  0,  // halo
//                          -7,7, -70,70, -700,700,  // core
//                          -8,8, -80,80, -800,800,  // core
//                          -9,9, -90,90, -900,900,  // core
//                          0,0,   0,10,    0,  0,  // halo
//                          0,0,   0,20,    0,  0}; // halo
//          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
//          break;
//        }
//      }
    }

    SECTION( "test_rank2_v2" )
    {
//      array::ArrayT<POD> arr(f.N,3,2);
//      array::ArrayView<POD,3> arrv = array::make_view<POD,3>(arr);
//      for( int p=0; p<f.N; ++p )
//      {
//        for( size_t i=0; i<3; ++i )
//        {
//          arrv(p,i,0) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 : -f.gidx[p]*std::pow(10,i) );
//          arrv(p,i,1) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 :  f.gidx[p]*std::pow(10,i) );
//        }
//      }

//      size_t strides[] = {2,2};
//      size_t shape[] = {3,1};
//      f.halo_exchange.execute(&arrv(0,0,1),strides,shape,2);

//      switch( parallel::mpi::comm().rank() )
//      {
//        case 0:
//        {
//          POD arr_c[] = {  0,9,   0,90,    0,900,  // halo
//                          -1,1, -10,10, -100,100,  // core
//                          -2,2, -20,20, -200,200,  // core
//                          -3,3, -30,30, -300,300,  // core
//                          0,4,   0,40,    0,400}; // halo
//          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
//          break;
//        }
//        case 1:
//        {
//          POD arr_c[] = {  0,3,   0,30,    0,300,  // halo
//                          -4,4, -40,40, -400,400,  // core
//                          -5,5, -50,50, -500,500,  // core
//                          -6,6, -60,60, -600,600,  // core
//                          0,7,   0,70,    0,700,  // halo
//                          0,8,   0,80,    0,800}; // halo
//          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
//          break;
//        }
//        case 2:
//        {
//          POD arr_c[] = {  0,5,   0,50,    0,500,  // halo
//                          0,6,   0,60,    0,600,  // halo
//                          -7,7, -70,70, -700,700,  // core
//                          -8,8, -80,80, -800,800,  // core
//                          -9,9, -90,90, -900,900,  // core
//                          0,1,   0,10,    0,100,  // halo
//                          0,2,   0,20,    0,200}; // halo
//          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
//          break;
//        }
//      }
    }

    SECTION( "test_rank0_ArrayView" )
    {
      array::ArraySpec spec(array::make_shape(1), array::make_strides(f.N));
      eckit::SharedPtr<array::Array> arr ( array::Array::wrap<POD>(f.gidx.data(), spec ) );
      array::ArrayView<POD,1> arrv = array::make_view<POD,1>(*arr);

      f.halo_exchange.execute(arrv);

      switch( parallel::mpi::comm().rank() )
      {
        case 0: { POD gidx_c[] = { 9, 1, 2, 3, 4};
          EXPECT(f.gidx == make_view(gidx_c,gidx_c+f.N)); break; }
        case 1: { POD gidx_c[] = { 3, 4, 5, 6, 7, 8};
          EXPECT(f.gidx == make_view(gidx_c,gidx_c+f.N)); break; }
        case 2: { POD gidx_c[] = { 5, 6, 7, 8, 9, 1, 2};
          EXPECT(f.gidx == make_view(gidx_c,gidx_c+f.N)); break; }
      }
    }


    SECTION( "test_rank1_ArrayView" )
    {
      array::ArrayT<POD> arr(f.N,2);
      array::ArrayView<POD,2> arrv = array::make_view<POD,2>(arr);
      for( int j=0; j<f.N; ++j ) {
        arrv(j,0) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j]*10 );
        arrv(j,1) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j]*100);
      }

      EXPECT( arr.stride(1) == 1 );
      // size_t strides[] = {1};
      // BOOST_CHECK_EQUAL_COLLECTIONS( arrv.strides()+1,arrv.strides()+2, strides,strides+1);

      f.halo_exchange.execute(arrv);

      switch( parallel::mpi::comm().rank() )
      {
        case 0: { POD arr_c[] = { 90,900, 10,100, 20,200, 30,300, 40,400 };
          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
        case 1: { POD arr_c[] = { 30,300, 40,400, 50,500, 60,600, 70,700, 80,800};
          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
        case 2: { POD arr_c[] = { 50,500, 60,600, 70,700, 80,800, 90,900, 10,100, 20,200};
          EXPECT(make_view(arrv.data(),arrv.data()+2*f.N) == make_view(arr_c,arr_c+2*f.N)); break; }
      }
    }

    SECTION( "test_rank2_ArrayView" )
    {
      array::ArrayT<POD> arr(f.N,3,2);
      array::ArrayView<POD,3> arrv = array::make_view<POD,3>(arr);
      for( int p=0; p<f.N; ++p )
      {
        for( size_t i=0; i<3; ++i )
        {
          arrv(p,i,0) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 : -f.gidx[p]*std::pow(10,i) );
          arrv(p,i,1) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 :  f.gidx[p]*std::pow(10,i) );
        }
      }

      EXPECT( arr.stride(1) == 2 );
      EXPECT( arr.stride(2) == 1 );
      // size_t strides[] = {2,1};
      // BOOST_CHECK_EQUAL_COLLECTIONS( arrv.strides()+1,arrv.strides()+3, strides,strides+2);

      f.halo_exchange.execute(arrv);

      switch( parallel::mpi::comm().rank() )
      {
        case 0:
        {
          POD arr_c[] = { -9,9, -90,90, -900,900,
                          -1,1, -10,10, -100,100,
                          -2,2, -20,20, -200,200,
                          -3,3, -30,30, -300,300,
                          -4,4, -40,40, -400,400};
          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
          break;
        }
        case 1:
        {
          POD arr_c[] = { -3,3, -30,30, -300,300,
                          -4,4, -40,40, -400,400,
                          -5,5, -50,50, -500,500,
                          -6,6, -60,60, -600,600,
                          -7,7, -70,70, -700,700,
                          -8,8, -80,80, -800,800};
          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
          break;
        }
        case 2:
        {
          POD arr_c[] = { -5,5, -50,50, -500,500,
                          -6,6, -60,60, -600,600,
                          -7,7, -70,70, -700,700,
                          -8,8, -80,80, -800,800,
                          -9,9, -90,90, -900,900,
                          -1,1, -10,10, -100,100,
                          -2,2, -20,20, -200,200};
          EXPECT(make_view(arrv.data(),arrv.data()+6*f.N) == make_view(arr_c,arr_c+6*f.N));
          break;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char **argv) {
    atlas::test::AtlasTestEnvironment env( argc, argv );
    return run_tests ( argc, argv, false );
}

