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
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using namespace std;
using namespace eckit;

//-----------------------------------------------------------------------------

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_field_pinned_mem") {
    atlas::util::Config config;
    config.set("host_memory_pinned", true);
    config.set("host_memory_mapped", false);
    auto field = Field("0", make_datatype<double>(), array::make_shape(10,4), config);

    auto view = array::make_view<double,2>(field);
    double* cpu_ptr = static_cast<double*>(view.data());
    view(3,2) = 1.;

#if ! ATLAS_HAVE_GRIDTOOLS_STORAGE
// TODO: gridtools storage does not implement view.index(...) at the moment

    cpu_ptr[view.index(3,2)] = 2.;
    EXPECT_EQ(view(3,2), 2.);

    field.updateDevice();
#pragma acc kernels present(cpu_ptr)
    {
        cpu_ptr[view.index(3,2)] = 3.;
    }

    field.updateHost();
    EXPECT_EQ(view(3,2), 3.);
#endif
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
