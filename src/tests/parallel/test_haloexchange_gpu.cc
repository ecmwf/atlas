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

CASE("test_haloexchange_gpu") {
  SETUP("Fixture") {
    Fixture f;

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
  }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char **argv) {
    atlas::test::AtlasTestEnvironment env( argc, argv );
    return run_tests ( argc, argv, false );
}

