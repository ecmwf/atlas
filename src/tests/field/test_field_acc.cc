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
#include "atlas/field/MultiField.h"

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


    auto dview = array::make_device_view<double,2>(field);
    double* dptr = dview.data();
#pragma acc parallel deviceptr(dptr)
    {
        dptr[dview.index(3,2)] = 4.;
    }
    field.updateHost();
    EXPECT_EQ( view(3,2), 4. );

#endif
}


CASE("test_fieldset_acc") {

    auto init_field_on_host = [](Field& field, double seed) {
        auto hview = array::make_host_view<double, 2>(field);
        for (int i=0; i<field.shape(0); ++i) {
            for (int j=0; j<field.shape(1); ++j) {
                hview(i,j) = seed + 100 + 10*i + j;
            }
        }
    };

    auto edit_field_on_device = [](Field& field) {
        int Ni = field.shape(0);
        int Nj = field.shape(1);

#if ATLAS_HAVE_GRIDTOOLS_STORAGE
        // TODO: gridtools storage does not implement view.index(...) at the moment
        {
            auto hview = array::make_host_view<double, 2>(field);
            for (int i=0; i<field.shape(0); ++i) {
                for (int j=0; j<field.shape(1); ++j) {
                    hview(i,j) += 10 + 1;
                }
            }
        }
#else
        auto dview = array::make_device_view<double, 2>(field);
        double* dptr = dview.data();
#pragma acc parallel deviceptr(dptr)
        {
            for (int i=0; i<Ni; ++i) {
                for (int j=0; j<Nj; ++j) {
                    dptr[dview.index(i,j)] += 10 + 1;
                }
            }
        }
#endif
    };
    auto verify_field_on_host = [](const Field& field, double seed) {
        auto hview = array::make_host_view<const double, 2>(field);
        for (int i=0; i<field.shape(0); ++i) {
            for (int j=0; j<field.shape(1); ++j) {
                EXPECT_EQ(hview(i,j), seed + 100 + 10 * (i+1) + (j+1));
            }
        }
    };

    auto create_fieldset = [] {
        const std::vector<int> vshape = std::vector<int>({4, -1, 3});
        const std::vector<std::string> var_names = {"temperature", "pressure", "density"};
        field::MultiField mfield(array::make_datatype<double>(), vshape, var_names);
        FieldSet fieldset;
        fieldset.add(mfield);
        fieldset.add(Field("contiguous", make_datatype<double>(), array::make_shape(4,3)));
        return fieldset;
    };

    SECTION("all fields in fieldset") {
        auto fieldset = create_fieldset();

        double seed = 0.;
        for( auto field : fieldset ) {
            init_field_on_host(field, seed);
            seed += 100.;
        }

        fieldset.allocateDevice();
        fieldset.updateDevice();
        for( auto field : fieldset ) {
            edit_field_on_device(field);
        }
        fieldset.updateHost();
        fieldset.deallocateDevice();

        seed = 0.;
        for( auto field : fieldset ) {
            verify_field_on_host(field, seed);
            seed += 100.;
        }
    }

    SECTION("field 'temperature' in fieldset") {
        auto fieldset = create_fieldset();
        auto field = fieldset.field("temperature");
        double seed = 0.;
        init_field_on_host(field, seed);
        fieldset.allocateDevice({"temperature"});
        fieldset.updateDevice({"temperature"});
        edit_field_on_device(field);
        fieldset.updateHost({"temperature"});
        fieldset.deallocateDevice({"temperature"});
        verify_field_on_host(field, seed);
    }

    SECTION("field [0] in fieldset") {
        auto fieldset = create_fieldset();
        auto field = fieldset.field(0);
        double seed = 0.;
        init_field_on_host(field, seed);
        fieldset.allocateDevice({0});
        fieldset.updateDevice({0});
        edit_field_on_device(field);
        fieldset.updateHost({0});
        fieldset.deallocateDevice({0});
        verify_field_on_host(field, seed);
    }
}
//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
