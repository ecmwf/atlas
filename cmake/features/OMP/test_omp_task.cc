/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


// This test mirrors functionality from src/atlas/parallel/omp/sort.h
// It appears that on some installations / compilers (such as apple-clang with llvm omp)
// the tasking does not seem to work properly.
// It all compiles but leads to runtime errors. This executable is
// added to allow compiler introspection for this feature.

#include <algorithm>
#include <vector>
#include <iostream>

#include <omp.h>


void recursive_task( int begin, int end ) {
    int size = end - begin;

    if ( size >= 2 ) {  // should be much larger in real case (e.g. 256)
        int mid = begin + size / 2;
        std::cout << begin << " - " << end  << " [" <<size << "]" << std::endl;
        {
#pragma omp task
            recursive_task( begin, mid );
#pragma omp task
            recursive_task( mid, end );
#pragma omp taskwait
        }
        // do work after coming back from nested recursion
    }
    else {
        // deepest level work
    }
}


void start_task( int size ) {
#pragma omp parallel
#pragma omp single
    recursive_task( 0, size );
}

int main() {
    start_task( 8 );
    return 0;
}
