/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cuda_runtime.h>
#include "tests/AtlasTestEnvironment.h"
#include "eckit/testing/Test.h"
#include "atlas/array.h"
#include "atlas/array/MakeView.h"
#include "atlas/runtime/Log.h"

using namespace atlas::array;
using namespace eckit::testing;

namespace atlas {
namespace test {

template<typename Value, int RANK>
__global__
void kernel_ex(ArrayView<Value, RANK> dv)
{
    dv(3, 3, 3) += dv.data_view().template length<0>() * dv.data_view().template length<1>() * dv.data_view().template length<2>();
}

CASE( "test_array" )
{
   constexpr unsigned int dx = 5;
   constexpr unsigned int dy = 6;
   constexpr unsigned int dz = 7;

   Array* ds = Array::create<double>(dx, dy, dz);
   ArrayView<double,3> hv = make_host_view<double, 3>(*ds);
   hv(3, 3, 3) = 4.5;

   ds->cloneToDevice();

   auto cv = make_device_view<double, 3>(*ds);
 
   kernel_ex<<<1,1>>>(cv);

   cudaDeviceSynchronize();

   ds->cloneFromDevice();
   ds->reactivateHostWriteViews();

   EXPECT( hv(3, 3, 3) == 4.5 + dx*dy*dz );

   delete ds;
}

}
}

int main(int argc, char **argv) {
    atlas::test::AtlasTestEnvironment env( argc, argv );
    return run_tests ( argc, argv, false );
}
