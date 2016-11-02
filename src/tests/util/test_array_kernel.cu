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

using namespace atlas::array;

namespace atlas {
namespace test {

template<typename Value, int RANK>
__global__
void kernel_ex(ArrayView<Value, RANK>& dv)
{
    dv(3) += 1;
}

BOOST_AUTO_TEST_CASE( test_array )
{
   Array* ds = Array::create<double>(4ul);
   auto hv = make_host_view<double, 1>(ds);
   hv(3) = 4.5;

   auto cv = make_device_view<double, 1>(ds);

   kernel_ex<<<1,1>>>(cv);

   cudaDeviceSynchronize();

   ds->reactivate_host_write_views();

   BOOST_CHECK_EQUAL( hv(3) , 5.5 );

   delete ds;
}

}
}
