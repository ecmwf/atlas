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

int test_throw() {
    std::cout << "--- " << __FUNCTION__ << std::endl;
    // For the dummy backend we expect to have an error thrown when used
    void* p;
    try {
        std::ignore = hicMalloc(&p, 8);
        std::ignore = hicFree(p);
    } catch(std::runtime_error& e) {
        std::string expected_message = "hicMalloc is using the dummy backend and should not be called";
        if (e.what() == expected_message) {
            std::cout << "--- " << __FUNCTION__ << " SUCCEEDED " << std::endl; 
            return 0; // success
        }
    }
    std::cout << "--- " << __FUNCTION__ << " SUCCEEDED " << std::endl; 
    return 1; // fail
}

std::vector<std::function<int()>> tests = { test_throw };

int main(int /*argc*/, char** /*argv*/) {
    int error = 0;
    for( auto& test: tests) {
        error += test();
    }
    return error;
}
