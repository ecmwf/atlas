/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "hic/hic.h"
#include "tests/AtlasTestEnvironment.h"
#include "atlas/array.h"
#include "atlas/array/MakeView.h"
#include "atlas/runtime/Log.h"

using namespace atlas::array;

namespace atlas {
namespace test {

#define REQUIRE_CUDA_SUCCESS(msg) \
do { \
  hicError_t err = hicPeekAtLastError(); \
  if (err != hicSuccess ) { \
    throw eckit::testing::TestException("REQUIRE_CUDA_SUCCESS ["+std::string(msg)+"] failed:\n"\
          + std::string(hicGetErrorString(err)) , Here()); \
  } \
} while(false)

template<typename Value, int RANK>
__global__
void kernel_ex(array::ArrayView<Value, RANK> dv)
{
    dv(3, 3, 3) += dv.shape(0) * dv.shape(1) * dv.shape(2);
}

template<typename Value, int RANK>
__global__
void loop_kernel_ex(array::ArrayView<Value, RANK> dv)
{
    for(int i=0; i < dv.shape(0); i++) {
      for(int j=0; j < dv.shape(1); j++) {
        for(int k=0; k < dv.shape(2); k++) {
          dv(i,j,k) += i*10+j*100+k*1000;
        }
      }
    }
}

CASE( "test_array" )
{
   constexpr unsigned int dx = 5;
   constexpr unsigned int dy = 6;
   constexpr unsigned int dz = 7;

   auto ds = std::unique_ptr<Array>(Array::create<double>(dx, dy, dz));
   auto hv = make_host_view<double, 3>(*ds);
   hv(3, 3, 3) = 4.5;

   ds->updateDevice();

   REQUIRE_CUDA_SUCCESS("updateDevice");

   auto cv = make_device_view<double, 3>(*ds);

   REQUIRE_CUDA_SUCCESS("make_device_view");

   kernel_ex<<<1,1>>>(cv);

   REQUIRE_CUDA_SUCCESS("kernel_ex");

   hicDeviceSynchronize();

   ds->updateHost();

   REQUIRE_CUDA_SUCCESS("updateHost");

   EXPECT( hv(3, 3, 3) == 4.5 + dx*dy*dz );
}

CASE( "test_array_loop" )
{
   constexpr unsigned int dx = 5;
   constexpr unsigned int dy = 6;
   constexpr unsigned int dz = 7;


   auto ds = std::unique_ptr<Array>(Array::create<double>(dx, dy, dz));
   array::ArrayView<double,3> hv = make_host_view<double, 3>(*ds);
   for(int i=0; i < dx; i++) {
     for(int j=0; j < dy; j++) {
       for(int k=0; k < dz; k++) {
         hv(i,j,k) = 0;
       }
     }
   }

   ds->updateDevice();

   REQUIRE_CUDA_SUCCESS("updateDevice");

   auto cv = make_device_view<double, 3>(*ds);

   REQUIRE_CUDA_SUCCESS("make_device_view");

   loop_kernel_ex<<<1,1>>>(cv);

   REQUIRE_CUDA_SUCCESS("loop_kernel_ex");

   hicDeviceSynchronize();

   ds->updateHost();

   REQUIRE_CUDA_SUCCESS("updateHost");

   for(int i=0; i < dx; i++) {
     for(int j=0; j < dy; j++) {
       for(int k=0; k < dz; k++) {
         EXPECT( hv(i,j,k) == i*10+j*100+k*1000 );
       }
     }
   }
}
}
}

int main(int argc, char **argv) {
    return atlas::test::run( argc, argv );
}
