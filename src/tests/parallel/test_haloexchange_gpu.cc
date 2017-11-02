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

template<typename DATA_TYPE, int Rank>
struct validate; 

template<typename DATA_TYPE>
struct validate<DATA_TYPE, 1> {

  static bool apply(array::ArrayView<DATA_TYPE,1>& arrv, DATA_TYPE arr_c[] ) {
    int strides[1];
    strides[0] = 1;
    for(size_t i = 0; i < arrv.template shape<0>(); ++i) {
std::cout << parallel::mpi::comm().rank() << "] VAL " << arrv(i) <<  " " << i << " " << arr_c[i*strides[0]] << std::endl;
        EXPECT(arrv(i) == arr_c[i*strides[0]]);
      
    }
  }
};

template<typename DATA_TYPE>
struct validate<DATA_TYPE, 2> {

  static bool apply(array::ArrayView<DATA_TYPE,2>& arrv, DATA_TYPE arr_c[] ) {
    int strides[2];
    strides[0] = 1;
    strides[1] = arrv.template shape<1>();
    for(size_t i = 0; i < arrv.template shape<0>(); ++i) {
      for(size_t j = 0; j < arrv.template shape<1>(); ++j) {
        EXPECT(arrv(i,j) == arr_c[i*strides[1] + j*strides[0]]);
      }
    }
  }
};


template<typename DATA_TYPE>
struct validate<DATA_TYPE, 3> {

  static bool apply(array::ArrayView<DATA_TYPE,3>& arrv, DATA_TYPE arr_c[] ) {
    int strides[3];
    strides[0] = 1;
    strides[1] = arrv.template shape<2>();
    strides[2] = arrv.template shape<1>()*strides[1];

    for(size_t i = 0; i < arrv.template shape<0>(); ++i) {
      for(size_t j = 0; j < arrv.template shape<1>(); ++j) {
        for(size_t k = 0; k < arrv.template shape<2>(); ++k) {
        EXPECT(arrv(i,j,k) == arr_c[i*strides[2] + j*strides[1] + k*strides[0]]);
        }
      }
    }
  }
};

CASE("test_haloexchange_gpu") {
  SETUP("Fixture") {
    Fixture f;

    SECTION( "test_rank0_arrview" )
    {
      array::ArrayT<POD> arr(f.N);
      array::ArrayView<POD,1> arrv = array::make_host_view<POD,1>(arr);
      for( int j=0; j<f.N; ++j ) {
        arrv(j) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j] );
      }

      arr.syncHostDevice();

      f.halo_exchange.template execute<POD,1>(arr, true);

      arr.syncHostDevice();

      switch( parallel::mpi::comm().rank() )
      {
        case 0: { POD arr_c[] = { 9, 1, 2, 3, 4};
          validate<POD,1>::apply(arrv, arr_c); break; }
        case 1: { POD arr_c[] = { 3, 4, 5, 6, 7, 8};
          validate<POD,1>::apply(arrv, arr_c); break; }
        case 2: { POD arr_c[] = { 5, 6, 7, 8, 9, 1, 2};
          validate<POD,1>::apply(arrv, arr_c); break; }
      }
      cudaDeviceSynchronize();
    }


    SECTION( "test_rank1" )
    {
      array::ArrayT<POD> arr(f.N,2);
      array::ArrayView<POD,2> arrv = array::make_host_view<POD,2>(arr);
      for( int j=0; j<f.N; ++j ) {
        arrv(j,0) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j]*10 );
        arrv(j,1) = (size_t(f.part[j]) != parallel::mpi::comm().rank() ? 0 : f.gidx[j]*100);
      }

      arr.syncHostDevice();

      f.halo_exchange.template execute<POD,2>(arr, true);

      arr.syncHostDevice();

      switch( parallel::mpi::comm().rank() )
      {
        case 0: { POD arr_c[] = { 90,900, 10,100, 20,200, 30,300, 40,400 };
          validate<POD,2>::apply(arrv, arr_c); break; }
        case 1: { POD arr_c[] = { 30,300, 40,400, 50,500, 60,600, 70,700, 80,800};
          validate<POD,2>::apply(arrv, arr_c); break; }
        case 2: { POD arr_c[] = { 50,500, 60,600, 70,700, 80,800, 90,900, 10,100, 20,200};
          validate<POD,2>::apply(arrv, arr_c); break; }
      }
    }

    SECTION( "test_rank2" )
    {
      array::ArrayT<POD> arr(f.N,3,2);
      array::ArrayView<POD,3> arrv = array::make_host_view<POD,3>(arr);
      for( int p=0; p<f.N; ++p )
      {
        for( size_t i=0; i<3; ++i )
        {
          arrv(p,i,0) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 : -f.gidx[p]*std::pow(10,i) );
          arrv(p,i,1) = (size_t(f.part[p]) != parallel::mpi::comm().rank() ? 0 :  f.gidx[p]*std::pow(10,i) );
        }
      }

      arr.syncHostDevice();

      f.halo_exchange.template execute<POD,3>(arr, true);

      arr.syncHostDevice();

      switch( parallel::mpi::comm().rank() )
      {
        case 0:
        {
          POD arr_c[] = { -9,9, -90,90, -900,900,
                          -1,1, -10,10, -100,100,
                          -2,2, -20,20, -200,200,
                          -3,3, -30,30, -300,300,
                          -4,4, -40,40, -400,400};
          validate<POD,3>::apply(arrv, arr_c); break; 
        }
        case 1:
        {
          POD arr_c[] = { -3,3, -30,30, -300,300,
                          -4,4, -40,40, -400,400,
                          -5,5, -50,50, -500,500,
                          -6,6, -60,60, -600,600,
                          -7,7, -70,70, -700,700,
                          -8,8, -80,80, -800,800};
          validate<POD,3>::apply(arrv, arr_c); break;
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
          validate<POD,3>::apply(arrv, arr_c); break; 
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

