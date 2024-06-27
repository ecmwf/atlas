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
#pragma acc parallel present(cpu_ptr)
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

void test_wrapping_a_slice(array::LocalView<double, 3> slice) {
  double* ptr = slice.data();
  array::ArrayShape shape(slice.shape(), slice.rank());
  array::ArrayStrides strides({slice.stride(0), slice.stride(1), slice.stride(2)});
  auto field = Field("name", ptr, array::ArraySpec(shape, strides));

  auto hview = array::make_host_view<double, 3>(field);
  for (idx_t i = 0; i < hview.shape(0); ++i) {
      for (idx_t j = 0; j < hview.shape(1); ++j) {
          for (idx_t k = 0; k < hview.shape(2); ++k) {
              hview(i,j,k) = 1000.*i + 100.*j + k;
          }
      }
  }

  field.updateDevice();
  auto dview = array::make_device_view<double,3>(field);
  double* dptr = dview.data();

#pragma acc kernels deviceptr(dptr)
  for (idx_t i=0; i < dview.shape(0); ++i) {
      for (idx_t j=0; j < dview.shape(1); ++j) {
          for (idx_t k=0; k < dview.shape(2); ++k) {
              dptr[dview.index(i,j,k)] *= -1;
          }
      }
  }

  // check host data before
  for (idx_t i=0; i < hview.shape(0); ++i) {
      for (idx_t j=0; j < hview.shape(1); ++j) {
          for (idx_t k=0; k < hview.shape(2); ++k) {
              EXPECT_EQ( hview(i,j,k), 1000.*i + 100.*j + k );
          }
      }
  }

  field.updateHost();
  field.deallocateDevice();

  // check host data after
  for (idx_t i = 0; i < hview.shape(0); ++i) {
      for (idx_t j = 0; j < hview.shape(1); ++j) {
          for (idx_t k = 0; k < hview.shape(2); ++k) {
              EXPECT_EQ( hview(i,j,k), -1000.*i - 100.*j - k );
          }
      }
  }
}

CASE("test_wrapping_discontiguous_data") {
  auto multifield = Field("name",make_datatype<double>(), array::make_shape(4,3,2,8));
  auto multiview = array::make_view<double,4>(multifield);
  multiview.assign(0.);

  auto all = array::Range::all();
  {
    auto slice = multiview.slice(2, all, all, all);
    test_wrapping_a_slice(slice);
  }
  {
    auto slice = multiview.slice(all, 0, all, all);
    test_wrapping_a_slice(slice);
  }
  {
    auto slice = multiview.slice(all, all, 1, all);
    test_wrapping_a_slice(slice);
  }
  {
    auto slice = multiview.slice(all, all, all, 6);
    test_wrapping_a_slice(slice);
  }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
