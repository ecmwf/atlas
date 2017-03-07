/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestArrayKernel
#include <cuda_runtime.h>
#include "ecbuild/boost_test_framework.h"
#include "atlas/array/Vector.h"
#include "atlas/array/gridtools/GPUClonable.h"


using namespace atlas::array;

namespace atlas {
namespace test {

struct int_gpu {
    int_gpu(int val) : val_(val), gpu_clone_(this) {}

    int_gpu* gpu_object_ptr() {return gpu_clone_.gpu_object_ptr();}

    void cloneToDevice(){ gpu_clone_.cloneToDevice();}

    int val_;
private:

    array::gridtools::GPUClonable<int_gpu> gpu_clone_;
};

__global__
void kernel_ex(VectorView<int_gpu*> list_ints)
{
    for(size_t i=0; i < list_ints.size(); ++i) {
        list_ints[i]->val_ += 5;
    }
}

BOOST_AUTO_TEST_CASE( test_vector_kernel )
{
    Vector<int_gpu*> list_ints(4);

    VectorView<int_gpu*> list_ints_h = make_host_vector_view(list_ints);

    list_ints_h[0] = new int_gpu(3);
    list_ints_h[1] = new int_gpu(4);
    list_ints_h[2] = new int_gpu(5);
    list_ints_h[3] = new int_gpu(6);

    list_ints.cloneToDevice();
    VectorView<int_gpu*> list_ints_d = make_device_vector_view(list_ints);

    kernel_ex<<<1,1>>>(list_ints_d);
}

}
}
