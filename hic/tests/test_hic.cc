/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#include <iostream>
#include <functional>
#include <vector>
#include <utility>

#include "hic/hic.h"

// -----------------------------------------------------------------------------

int test_hicMalloc() {
    std::cout << "--- " << __FUNCTION__ << std::endl;
    void* p;
    HIC_CALL( hicMalloc(&p, 8) );
    HIC_CALL( hicFree(p) );
    std::cout << "--- " << __FUNCTION__ << " SUCCEEDED " << std::endl; 
    return 0; // success 
}

// -----------------------------------------------------------------------------

int test_hicMallocManaged() {
    std::cout << "--- " << __FUNCTION__ << std::endl;
    void* p;
    HIC_CALL( hicMallocManaged(&p, 8) );
    HIC_CALL( hicFree(p) );
    std::cout << "--- " << __FUNCTION__ << " SUCCEEDED " << std::endl; 
    return 0; // success 
}

// -----------------------------------------------------------------------------

std::vector<std::function<int()>> tests = {
    test_hicMalloc,
    test_hicMallocManaged,
};

int main(int argc, char* argv[]) {
    int num_devices = 0;
    hicGetDeviceCount(&num_devices);
    if( num_devices == 0 ) {
        std::ignore = hicGetLastError();
        std::cout << "TEST IGNORED, hicGetDeviceCount -> 0" << std::endl; 
        return 0;
    }
    std::cout << "hicGetDeviceCount -> " << num_devices << std::endl; 
    int error = 0;
    for( auto& test: tests) {
        try {
            error += test();
        }
        catch( std::exception& e ) {
            error += 1;
            std::cout << e.what() << std::endl;
        }
    }
    return error;
}
