/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

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
#ifndef _OPENACC
#error No OpenACC
#endif

    int* dev;
    cudaMalloc(&dev, sizeof(int));
    int* c_ptr = new int();
    *c_ptr = 5;

    acc_map_data(c_ptr, dev, sizeof(int));
#pragma acc kernels present(c_ptr)
    *c_ptr = 2.;

    int* h = new int();
    cudaMemcpy(h, dev, sizeof(int), cudaMemcpyDeviceToHost);
    EXPECT_EQ( *h, 2 );
}


CASE("test_field_acc") {
    auto field = Field("0", make_datatype<double>(), array::make_shape(10,4));

    auto view = array::make_view<double,2>(field);
    double* cpu_ptr = static_cast<double*>(view.data());
    view(3,2) = 1.;
    cpu_ptr[view.index(3,2)] = 2.;

    EXPECT_EQ( view(3,2), 2. );

    field.updateDevice();

#pragma acc kernels present(cpu_ptr)
    {
        cpu_ptr[view.index(3,2)] = 3.;
    }

    field.updateHost();

    EXPECT_EQ( view(3,2), 3. );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
