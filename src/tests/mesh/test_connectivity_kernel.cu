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
void kernel_ex(BlockConnectivity* conn, bool* result)
{

    *result &= (conn->rows() == 2);
    *result &= (conn->cols() == 5);

    *result &= ((*conn)(0,2) == 9 + FROM_FORTRAN+FORTRAN_BASE);
    *result &= ((*conn)(0,4) == 356 + FROM_FORTRAN+FORTRAN_BASE);
    *result &= ((*conn)(1,1) == 3 + FROM_FORTRAN+FORTRAN_BASE);

}

BOOST_AUTO_TEST_CASE( test_connectivity )
{
    BlockConnectivity* conn;
    cudaMallocManaged(&conn, sizeof(BlockConnectivity));

    bool* result;
    cudaMallocManaged(&result, sizeof(bool));

    *result = true;

    idx_t vals2[12] = {2,3,9,34,356,86,3,24,84,45,2,2};

    conn->add(2,5, vals2);


//    Array* ds = Array::create<double>(4ul, 4ul, 4ul);
//   ArrayView<double,3> hv = make_host_view<double, 3>(*ds);
//   hv(3, 3, 3) = 4.5;


//   ds->clone_to_device();


//   auto cv = make_device_view<double, 3>(*ds);

   kernel_ex<<<1,1>>>(conn, result);

   cudaDeviceSynchronize();

   BOOST_CHECK_EQUAL( *result , true );


//   ds->clone_from_device();
//   ds->reactivate_host_write_views();

//   BOOST_CHECK_EQUAL( hv(3, 3, 3) , 5.5 );

//   delete ds;
}

}
}
