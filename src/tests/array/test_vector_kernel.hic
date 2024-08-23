/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "hic/hic.h"

#include "tests/AtlasTestEnvironment.h"

#include "atlas/array/Vector.h"
#include "atlas/array.h"
#include "atlas/util/GPUClonable.h"
#include "atlas/array/MakeView.h"
#include "atlas/runtime/Log.h"

using namespace atlas::array;

namespace atlas {
namespace test {

struct int_gpu {
    int_gpu() = delete;
    int_gpu(int val) : val_(val), gpu_clone_(this) {}

    int_gpu* gpu_object_ptr() {return gpu_clone_.gpu_object_ptr();}

    void updateDevice(){ gpu_clone_.updateDevice();}
    void updateHost(){ gpu_clone_.updateHost();}

    int val_;

private:
    util::GPUClonable<int_gpu> gpu_clone_;
};

__global__
void kernel_ex(VectorView<int_gpu*>* list_ints )
{
    for(size_t i=0; i < list_ints->size() ; ++i) {
        (*list_ints)[i]->val_ += 5;
    }
}

CASE( "test_resize" )
{
    Vector<int_gpu*> list_ints(2);

    VectorView<int_gpu*> list_ints_h = make_host_vector_view(list_ints);
    list_ints_h[0] = new int_gpu(3);
    list_ints_h[1] = new int_gpu(4);

    EXPECT( list_ints_h[0]->val_ == 3 );
    EXPECT( list_ints_h[1]->val_ == 4 );

    list_ints.resize(6);
    EXPECT( list_ints_h.is_valid(list_ints) == false );

    VectorView<int_gpu*> list_ints_h2 = make_host_vector_view(list_ints);


    EXPECT( list_ints_h2[0]->val_ == 3 );
    EXPECT( list_ints_h2[1]->val_ == 4 );
    EXPECT( list_ints_h2.size()   == 6 );
}



CASE( "test_vector_kernel" )
{
    Vector<int_gpu*> list_ints(4);

    VectorView<int_gpu*> list_ints_h = make_host_vector_view(list_ints);

    list_ints_h[0] = new int_gpu(3);
    list_ints_h[1] = new int_gpu(4);
    list_ints_h[2] = new int_gpu(5);
    list_ints_h[3] = new int_gpu(6);
    list_ints.updateDevice();

    VectorView<int_gpu*> list_ints_d = make_device_vector_view(list_ints);

    VectorView<int_gpu*>* list_ints_dp;
    hicMalloc((void**)(&list_ints_dp), sizeof(VectorView<int_gpu*>));

    hicMemcpy(list_ints_dp, &list_ints_d, sizeof(VectorView<int_gpu*>), hicMemcpyHostToDevice);

    kernel_ex<<<1,1>>>(list_ints_dp);

    if( hicPeekAtLastError() != hicSuccess) std::cout << "ERROR " << std::endl;
    list_ints.updateHost();

    EXPECT( list_ints_h[0]->val_ == 8 );

    EXPECT_THROWS_AS( list_ints.resize(6), eckit::AssertionFailed );

}

}
}

int main(int argc, char **argv) {
    return atlas::test::run( argc, argv );
}
