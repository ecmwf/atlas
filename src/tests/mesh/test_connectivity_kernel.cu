/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestConnectivityKernel
#include <cuda_runtime.h>
#include "ecbuild/boost_test_framework.h"
#include "atlas/mesh/Connectivity.h"

using namespace atlas::mesh;

namespace atlas {
namespace test {

#ifdef ATLAS_HAVE_FORTRAN
#define IN_FORTRAN -1
#else
#define IN_FORTRAN
#endif


__global__
void kernel_block(BlockConnectivityImpl* conn_, bool* result)
{
    BlockConnectivityImpl& conn = *conn_;

    *result &= (conn.rows() == 2);
    *result &= (conn.cols() == 5);

    *result &= ((conn)(0,2) == 9 );
    *result &= ((conn)(0,4) == 356 );
    *result &= ((conn)(1,1) == 3 );
}

__global__
void kernel_irr(IrregularConnectivityImpl* conn_, bool* result)
{

    IrregularConnectivityImpl& conn = *conn_;

    *result = true;

    *result &= (conn.rows()== 2);
    *result &= (conn.cols(0) == 3);
    *result &= (conn.cols(1) == 3);
    *result &= (conn.mincols() == 3);
    *result &= (conn.maxcols() == 3);

    *result &= (conn(0,0) == 1 IN_FORTRAN);
    *result &= (conn(0,1) == 3 IN_FORTRAN);
    *result &= (conn(0,2) == 4 IN_FORTRAN);

}

__global__
void kernel_multiblock(MultiBlockConnectivityImpl* conn_, bool* result)
{

    MultiBlockConnectivityImpl& conn = *conn_;

    *result = true;

    *result &= (conn.blocks()== 1);
    *result &= (conn.rows()== 2);
    *result &= (conn.cols(0) == 3);
    *result &= (conn.cols(1) == 3);
    *result &= (conn.mincols() == 3);
    *result &= (conn.maxcols() == 3);

    *result &= (conn(0,0) == 1 IN_FORTRAN);
    *result &= (conn(0,1) == 3 IN_FORTRAN);
    *result &= (conn(0,2) == 4 IN_FORTRAN);

    BlockConnectivity& ablock = conn.block(0);
    *result &= (ablock(0,0) == 1 IN_FORTRAN);
    *result &= (ablock(0,1) == 3 IN_FORTRAN);
    *result &= (ablock(0,2) == 4 IN_FORTRAN);
}



BOOST_AUTO_TEST_CASE( test_block_connectivity )
{
    BlockConnectivity conn;

    bool* result;
    cudaMallocManaged(&result, sizeof(bool));

    *result = true;

    idx_t vals2[12] = {2,3,9,34,356,86,3,24,84,45,2,2};

    conn.add(2,5, vals2);
    conn.cloneToDevice();
    BOOST_CHECK( !conn.deviceNeedsUpdate() );

    kernel_block<<<1,1>>>(conn.gpu_object_ptr(), result);

    cudaDeviceSynchronize();

    BOOST_CHECK_EQUAL( *result , true );

    // copy back, although not strickly needed since the gpu copy does not modify values, 
    // but for the sake of testing it

    conn.cloneFromDevice();
    BOOST_CHECK_EQUAL((conn)(0,4), 356 );
}

BOOST_AUTO_TEST_CASE( test_irregular_connectivity )
{
    IrregularConnectivity conn("mesh");
    BOOST_CHECK_EQUAL(conn.rows(),0);
    BOOST_CHECK_EQUAL(conn.maxcols(),0);

    constexpr idx_t vals[6] = {1,3,4,3,7,8};
    bool from_fortran = true;
    conn.add(2, 3, vals, from_fortran);

    bool* result;
    cudaMallocManaged(&result, sizeof(bool));
    *result = true;

    conn.cloneToDevice();
    BOOST_CHECK( !conn.deviceNeedsUpdate() );

    kernel_irr<<<1,1>>>(conn.gpu_object_ptr(), result);

    cudaDeviceSynchronize();

    BOOST_CHECK_EQUAL( *result , true );

    // copy back, although not strickly needed since the gpu copy does not modify values, 
    // but for the sake of testing it
    conn.cloneFromDevice();
    BOOST_CHECK_EQUAL(conn(0,1), 3 IN_FORTRAN);

}

BOOST_AUTO_TEST_CASE( test_multiblock_connectivity )
{
    
    MultiBlockConnectivity conn("mesh");
    BOOST_CHECK_EQUAL(conn.rows(),0);
    BOOST_CHECK_EQUAL(conn.maxcols(),0);

    constexpr idx_t vals[6] = {1,3,4,3,7,8};
    bool from_fortran = true;
    conn.add(2, 3, vals, from_fortran);
    
    BOOST_CHECK(conn.block(0)(0,0) == 1);
    bool* result;
    cudaMallocManaged(&result, sizeof(bool));
    *result = true;

    conn.cloneToDevice();
    BOOST_CHECK( !conn.deviceNeedsUpdate() );

    kernel_multiblock<<<1,1>>>(conn.gpu_object_ptr(), result);

    cudaDeviceSynchronize();

    BOOST_CHECK_EQUAL( *result , true );

    // copy back, although not strickly needed since the gpu copy does not modify values, 
    // but for the sake of testing it
    conn.cloneFromDevice();
    BOOST_CHECK(conn.block(0)(0,0) == 1);

}


}
}
