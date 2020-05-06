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
// It appears that on some installations / compilers the "untied if" in the
// openmp pragma compiles but leads to runtime errors. This executable is
// added to allow compiler introspection for this feature.

#include <algorithm>
#include <vector>
#include <omp.h>

template <typename RandomAccessIterator>
void merge_sort_recursive( const RandomAccessIterator& iterator, size_t begin, size_t end ) {
    auto size = end - begin;
    if ( size >= 2 ) { // should be much larger in real case (e.g. 256)
        auto mid = begin + size / 2;
        {
#pragma omp task shared( iterator ) untied if ( size >= ( 1 << 15 ) )
            merge_sort_recursive( iterator, begin, mid );
#pragma omp task shared( iterator ) untied if ( size >= ( 1 << 15 ) )
            merge_sort_recursive( iterator, mid, end );
#pragma omp taskwait
        }
        std::inplace_merge( iterator + begin, iterator + mid, iterator + end );
    }
    else {
        std::sort( iterator + begin, iterator + end );
    }
}

template <typename RandomAccessIterator>
void omp_sort( RandomAccessIterator first, RandomAccessIterator last ) {
#pragma omp parallel
#pragma omp single
        merge_sort_recursive( first, 0, std::distance( first, last ) );
}

int main() {
    auto integers = std::vector<int>(8);
    omp_sort( integers.begin(), integers.end() );
    return 0;
}
