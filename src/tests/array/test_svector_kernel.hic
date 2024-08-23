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

#include "atlas/library/config.h"
#include "tests/AtlasTestEnvironment.h"
#include "atlas/array/SVector.h"

using namespace atlas::array;

namespace atlas {
namespace test {

__global__
void kernel_exe(int* list_ints_ptr, size_t size, int offset, bool* result )
{
    SVector<int> list_ints(list_ints_ptr, size);

    *result = *result && (list_ints[offset] == 3);
    *result = *result && (list_ints[offset+1] == 4);


    list_ints[offset]++;
    list_ints[offset+1]++;
}

CASE( "test_svector" )
{
    SVector<int> list_ints(2);

    list_ints[0] = 3;
    list_ints[1] = 4;

    EXPECT( list_ints[0] == 3 );
    EXPECT( list_ints[1] == 4 );

    EXPECT( list_ints.size() == 2);

    bool *result;
    hicError_t err = hicMallocManaged(&result, sizeof(bool));

    if(err != hicSuccess)
        throw_AssertionFailed("failed to allocate GPU memory");

    *result=true;
    kernel_exe<<<1,1>>>(list_ints.data(), list_ints.size(), 0, result);
    hicDeviceSynchronize();

    err = hicGetLastError();
    if(err != hicSuccess)
        throw_AssertionFailed("failed to execute kernel");

    EXPECT( *result );
    EXPECT( list_ints[0] == 4);
    EXPECT( list_ints[1] == 5);

}

CASE( "test_svector_resize" )
{
    SVector<int> list_ints(2);

    list_ints[0] = 3;
    list_ints[1] = 4;

    EXPECT( list_ints[0] == 3 );
    EXPECT( list_ints[1] == 4 );

    EXPECT( list_ints.size() == 2);

    list_ints.resize(5);

    EXPECT( list_ints.size() == 5);

    bool *result;
    hicError_t err = hicMallocManaged(&result, sizeof(bool));

    if(err != hicSuccess)
        throw_AssertionFailed("failed to allocate GPU memory");

    *result=true;

    list_ints[3] = 3;
    list_ints[4] = 4;

    kernel_exe<<<1,1>>>(list_ints.data(), list_ints.size(), 3, result);
    hicDeviceSynchronize();

    err = hicGetLastError();
    if(err != hicSuccess)
        throw_AssertionFailed("failed to execute kernel");

    EXPECT( *result );
    EXPECT( list_ints[3] == 4);
    EXPECT( list_ints[4] == 5);

}

}
}

int main(int argc, char **argv) {
    return atlas::test::run( argc, argv );
}
