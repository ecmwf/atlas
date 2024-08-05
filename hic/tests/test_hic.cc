#include <iostream>
#include <functional>
#include <vector>

#include "hic/hic.h"

int test_hic() {
    std::cout << "--- " << __FUNCTION__ << std::endl;
    void* p;
    HIC_CALL( hicMalloc(&p, 8) );
    HIC_CALL( hicFree(p) );
    std::cout << "--- " << __FUNCTION__ << " SUCCEEDED " << std::endl; 
    return 0; // success 
}

std::vector<std::function<int()>> tests = { test_hic };

int main(int argc, char* argv[]) {
    int error = 0;
    for( auto& test: tests) {
        error += test();
    }
    return error;
}
