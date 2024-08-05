#include <iostream>
#include <functional>

#include "hic/hic.h"

int test_throw() {
    std::cout << "--- " << __FUNCTION__ << std::endl;
    // For the dummy backend we expect to have an error thrown when used
    void* p;
    try {
        HIC_CALL( hicMalloc(&p, 8) );
        HIC_CALL( hicFree(p) );
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

int main(int argc, char* argv[]) {
    int error = 0;
    for( auto& test: tests) {
        error += test();
    }
    return error;
}