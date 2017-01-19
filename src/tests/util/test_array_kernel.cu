/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestArrayKernel
#include <cuda_runtime.h>
#include "ecbuild/boost_test_framework.h"
#include "atlas/array/Array.h"
#include "atlas/array/MakeView.h"
#include "atlas/runtime/Log.h"

using namespace atlas::array;

namespace atlas {
namespace test {

template<typename Value, int RANK>
__global__
void kernel_ex(ArrayView<Value, RANK> dv)
{
    dv(3, 3, 3) += 1;
}

BOOST_AUTO_TEST_CASE( test_array )
{
   ATLAS_DEBUG_HERE();
   Array* ds = Array::create<double>(4ul, 4ul, 4ul);
   ArrayView<double,3> hv = make_host_view<double, 3>(*ds);
   hv(3, 3, 3) = 4.5;

   ATLAS_DEBUG_HERE();

   ds->cloneToDevice();

   ATLAS_DEBUG_HERE();

   auto cv = make_device_view<double, 3>(*ds);

   ATLAS_DEBUG_HERE();

   kernel_ex<<<1,1>>>(cv);

   ATLAS_DEBUG_HERE();

   cudaDeviceSynchronize();

   ATLAS_DEBUG_HERE();

   ds->cloneFromDevice();
   ds->reactivateHostWriteViews();

   BOOST_CHECK_EQUAL( hv(3, 3, 3) , 5.5 );

   delete ds;
}

}
}
