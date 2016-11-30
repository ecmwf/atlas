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
#define FORTRAN_BASE 1
#define INDEX_REF Index
#define FROM_FORTRAN -1
#define TO_FORTRAN   +1
#else
#define FORTRAN_BASE 0
#define INDEX_REF *
#define FROM_FORTRAN
#define TO_FORTRAN
#endif


__global__
void kernel_block(BlockConnectivity conn, bool* result)
{

    *result &= (conn.rows() == 2);
    *result &= (conn.cols() == 5);

    *result &= ((conn)(0,2) == 9 + FROM_FORTRAN+FORTRAN_BASE);
    *result &= ((conn)(0,4) == 356 + FROM_FORTRAN+FORTRAN_BASE);
    *result &= ((conn)(1,1) == 3 + FROM_FORTRAN+FORTRAN_BASE);
}

__global__
void kernel_irr(IrregularConnectivity conn, bool* result)
{


    *result = true;
    
    *result &= (conn.rows()== 2);
    *result &= (conn.cols(0) == 3);
    *result &= (conn.cols(1) == 3);
    *result &= (conn.mincols() == 3);
    *result &= (conn.maxcols() == 3);
printf(" %d %d %d ", conn(0,0), conn(0,1), conn(0,2));
    *result &= (conn(0,0),1 + FROM_FORTRAN);
    *result &= (conn(0,1),3 + FROM_FORTRAN);
    *result &= (conn(0,2),4 + FROM_FORTRAN);

}


BOOST_AUTO_TEST_CASE( test_block_connectivity )
{
    BlockConnectivity conn;

    bool* result;
    cudaMallocManaged(&result, sizeof(bool));

    *result = true;

    idx_t vals2[12] = {2,3,9,34,356,86,3,24,84,45,2,2};

    conn.add(2,5, vals2);
    conn.clone_to_device();
    kernel_block<<<1,1>>>(conn, result);

   cudaDeviceSynchronize();
 
   BOOST_CHECK_EQUAL( *result , true );

}

BOOST_AUTO_TEST_CASE( test_irregular_connectivity )
{
    IrregularConnectivity conn("mesh");
    BOOST_CHECK_EQUAL(conn.rows(),0);
    BOOST_CHECK_EQUAL(conn.maxcols(),0);

    constexpr idx_t vals[6] = {1,3,4,3,7,8};
    conn.add(2, 3, vals, true);

    bool* result;
    cudaMallocManaged(&result, sizeof(bool));

    *result = true;
    conn.clone_to_device();
    kernel_irr<<<1,1>>>(conn, result);

    cudaDeviceSynchronize();

    BOOST_CHECK_EQUAL( *result , true );

}

}
}
