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

#include "hic/hic.h"
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
    HIC_CALL(hicMalloc(&d_ptr, sizeof(int)));
    acc_map_data(c_ptr, d_ptr, sizeof(int));

    HIC_CALL(hicMemcpy(d_ptr, c_ptr, sizeof(int), hicMemcpyHostToDevice));

#pragma acc kernels present(c_ptr)
    {
        *c_ptr -= 3.;
    }

    EXPECT_EQ( *c_ptr, 5. );

    HIC_CALL(hicMemcpy(c_ptr, d_ptr, sizeof(int), hicMemcpyDeviceToHost));
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
#endif
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
