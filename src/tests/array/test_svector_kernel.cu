/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestSVectorKernel
#include <cuda_runtime.h>
#include "ecbuild/boost_test_framework.h"
#include "atlas/array/SVector.h"

using namespace atlas::array;

namespace atlas {
namespace test {

__global__
void kernel_exe(array::SVector<int> list_ints, int offset, bool* result )
{
    *result = *result && (list_ints[offset] == 3);
    *result = *result && (list_ints[offset+1] == 4);


    list_ints[offset]++;
    list_ints[offset+1]++;
}

BOOST_AUTO_TEST_CASE( test_svector )
{
    SVector<int> list_ints(2);

    list_ints[0] = 3;
    list_ints[1] = 4;

    BOOST_CHECK_EQUAL( list_ints[0] , 3 );
    BOOST_CHECK_EQUAL( list_ints[1] , 4 );

    BOOST_CHECK_EQUAL( list_ints.size(), 2);

    bool *result;
    cudaError_t err = cudaMallocManaged(&result, sizeof(bool));

    if(err != cudaSuccess)
        throw eckit::AssertionFailed("failed to allocate GPU memory");

    *result=true;
    kernel_exe<<<1,1>>>(list_ints, 0, result);
    cudaDeviceSynchronize();
    
    err = cudaGetLastError();
    if(err != cudaSuccess)
        throw eckit::AssertionFailed("failed to execute kernel");

    BOOST_CHECK( *result );
    BOOST_CHECK_EQUAL( list_ints[0], 4);
    BOOST_CHECK_EQUAL( list_ints[1], 5);

}

BOOST_AUTO_TEST_CASE( test_svector_resize )
{
    SVector<int> list_ints(2);

    list_ints[0] = 3;
    list_ints[1] = 4;

    BOOST_CHECK_EQUAL( list_ints[0] , 3 );
    BOOST_CHECK_EQUAL( list_ints[1] , 4 );

    BOOST_CHECK_EQUAL( list_ints.size(), 2);

    list_ints.resize(5);

    BOOST_CHECK_EQUAL( list_ints.size(), 5);

    bool *result;
    cudaError_t err = cudaMallocManaged(&result, sizeof(bool));

    if(err != cudaSuccess)
        throw eckit::AssertionFailed("failed to allocate GPU memory");

    *result=true;

    list_ints[3] = 3;
    list_ints[4] = 4;

    kernel_exe<<<1,1>>>(list_ints, 3, result);
    cudaDeviceSynchronize();

    err = cudaGetLastError();
    if(err != cudaSuccess)
        throw eckit::AssertionFailed("failed to execute kernel");

    BOOST_CHECK( *result );
    BOOST_CHECK_EQUAL( list_ints[3], 4);
    BOOST_CHECK_EQUAL( list_ints[4], 5);

}

}
}
