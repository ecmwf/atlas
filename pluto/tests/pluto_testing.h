
/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstddef>
#include <iostream>

int FAILED_EXPECTS = 0;
#define EXPECT(COND) \
    do { \
        if (! COND) { \
            std::cerr << "EXPECT failed: " << #COND << std::endl; \
            FAILED_EXPECTS += 1; \
        } \
    } while(false)

inline int pluto_testing_return() {
    if (FAILED_EXPECTS) {
        std::cerr << "TEST FAILED!  Failures: " << FAILED_EXPECTS << std::endl;
        return 1;
    }
    else {
        std::cerr << "TEST SUCCESS!" << std::endl;
        return 0;
    }
}
