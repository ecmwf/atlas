/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef _OPENACC
#error This file needs to be compiled with OpenACC support
#endif

#include <cuda_runtime.h>
#include <openacc.h>

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

#include "tests/AtlasTestEnvironment.h"

using namespace std;
using namespace eckit;

//-----------------------------------------------------------------------------

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_acc") {
    int* c_ptr = new int();
    *c_ptr = 5;

    int* d_ptr;
    cudaMalloc(&d_ptr, sizeof(int));
    acc_map_data(c_ptr, d_ptr, sizeof(int));

    cudaMemcpy(d_ptr, c_ptr, sizeof(int), cudaMemcpyHostToDevice);

#pragma acc kernels present(c_ptr)
    {
        *c_ptr -= 3.;
    }

    EXPECT_EQ( *c_ptr, 5. );

    cudaMemcpy(c_ptr, d_ptr, sizeof(int), cudaMemcpyDeviceToHost);
    EXPECT_EQ( *c_ptr, 2. );
}


CASE("test_field_acc") {
    auto field = Field("0", make_datatype<double>(), array::make_shape(10,4));

    auto view = array::make_view<double,2>(field);
    double* cpu_ptr = static_cast<double*>(view.data());
    view(3,2) = 1.;

#if ! ATLAS_HAVE_GRIDTOOLS_STORAGE
// TODO: gridtools storage does not implement view.index(...) at the moment

    cpu_ptr[view.index(3,2)] = 2.;
    EXPECT_EQ( view(3,2), 2. );

    field.updateDevice();
#pragma acc kernels present(cpu_ptr)
    {
        cpu_ptr[view.index(3,2)] = 3.;
    }
    field.updateHost();
    EXPECT_EQ( view(3,2), 3. );


    //auto dview = array::make_device_view<double,2>(field);
    auto dview = field.array().device_data<double>();
#pragma acc kernels deviceptr(dview)
    {
        //dview(3,2) = 4.;
        dview[view.index(3,2)] = 4.;
    }
    field.updateHost();
    EXPECT_EQ( view(3,2), 4. );

#endif
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
