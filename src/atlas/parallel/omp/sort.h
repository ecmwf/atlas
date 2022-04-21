/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <algorithm>
#include <functional>
#include <iterator>

#include "atlas/library/config.h"
#include "atlas/parallel/omp/omp.h"

#if ATLAS_HAVE_OMP && ATLAS_OMP_TASK_SUPPORTED
#include <omp.h>
#define ATLAS_HAVE_OMP_SORTING 1
#else
#define ATLAS_HAVE_OMP_SORTING 0
#endif

// Bug in Cray 8.5 or below results in segmentation fault in atlas_test_omp_sort
#if ATLAS_HAVE_OMP_SORTING && defined(_CRAYC)
#if _RELEASE <= 8 && _RELEASE_MINOR < 6
#undef ATLAS_HAVE_OMP_SORTING
#define ATLAS_HAVE_OMP_SORTING 0
#endif
#endif

namespace atlas {
namespace omp {

/**
 * sort
 * ====
 *
 * 1)  template <typename RandomAccessIterator>
 *       void sort ( RandomAccessIterator first, RandomAccessIterator last );
 *
 * 2)  template <typename RandomAccessIterator, typename Compare>
 *       void sort ( RandomAccessIterator first, RandomAccessIterator last, Compare comp );
 *
 * Sort elements in range
 * Sorts the elements in the range [first,last) into ascending order.
 *
 * The elements are compared using operator< for the first version, and comp for the second.
 *
 * Equivalent elements are not guaranteed to keep their original relative order (see stable_sort).
 *
 * Parameters
 * ----------
 * first, last
 *     Random-access iterators to the initial and final positions of the sequence to be sorted. The range used is [first,last),
 *     which contains all the elements between first and last, including the element pointed by first but not the element pointed by last.
 *     RandomAccessIterator shall point to a type for which swap is properly defined and which is both move-constructible and move-assignable.
 * comp
 *     Binary function that accepts two elements in the range as arguments, and returns a value convertible to bool.
 *     The value returned indicates whether the element passed as first argument is considered to go before the second in the specific strict weak ordering it defines.
 *     The function shall not modify any of its arguments.
 *     This can either be a function pointer or a function object.
 *
 *
 *
 * merge_blocks
 * ============
 *
 * 1)  template<typename RandomAccessIterator, typename RandomAccessIterator2>
 *       void merge_blocks( RandomAccessIterator first, RandomAccessIterator last,
 *                          RandomAccessIterator2 blocks_size_first, RandomAccessIterator2 blocks_size_last );
 *
 *
 * 1)  template<typename RandomAccessIterator, typename RandomAccessIterator2, typename Compare >
 *       void merge_blocks( RandomAccessIterator first, RandomAccessIterator last,
 *                          RandomAccessIterator2 blocks_size_first, RandomAccessIterator2 blocks_size_last,
 *                          Compare compare );
 *
 * Sort elements in range [first + *blocks_first, first + *blocks_last) using a merge sort
 * where each block in range [blocks_first,blocks_last) is already sorted.
 *
 * Parameters
 * ----------
 * first, last
 *     Random-access iterators for bounding the sequence to be sorted
 * blocks_begin, blocks_end
 *     Random-access iterators that define offsets from paramter "first" of blocks that are already sorted
 */


namespace detail {

#if ATLAS_HAVE_OMP_SORTING
template <typename RandomAccessIterator, typename Compare>
void merge_sort_recursive(const RandomAccessIterator& iterator, size_t begin, size_t end, Compare compare) {
    auto size = end - begin;
    if (size >= 256) {
        auto mid = begin + size / 2;
        //#pragma omp taskgroup
        //   --> it would be preferred to use taskgroup and taskyield instead of taskwait,
        //       but this leads to segfaults on Cray (cce/8.5.8)
        {
#if ATLAS_OMP_TASK_UNTIED_SUPPORTED
#pragma omp task shared(iterator) untied if (size >= (1 << 15))
#else
#pragma omp task shared(iterator)
#endif
            merge_sort_recursive(iterator, begin, mid, compare);
#if ATLAS_OMP_TASK_UNTIED_SUPPORTED
#pragma omp task shared(iterator) untied if (size >= (1 << 15))
#else
#pragma omp task shared(iterator)
#endif
            merge_sort_recursive(iterator, mid, end, compare);
//#pragma omp taskyield
#pragma omp taskwait
        }
        std::inplace_merge(iterator + begin, iterator + mid, iterator + end, compare);
    }
    else {
        std::sort(iterator + begin, iterator + end, compare);
    }
}
#endif

#if ATLAS_HAVE_OMP_SORTING
template <typename RandomAccessIterator, typename Indexable, typename Compare>
void merge_blocks_recursive(const RandomAccessIterator& iterator, const Indexable& blocks, size_t blocks_begin,
                            size_t blocks_end, Compare compare) {
    if (blocks_end <= blocks_begin + 1) {
        // recursion done, go back out
        return;
    }
    size_t blocks_mid = (blocks_begin + blocks_end) / 2;
    //#pragma omp taskgroup
    //   --> it would be preferred to use taskgroup and taskyield instead of taskwait,
    //       but this leads to segfaults on Cray (cce/8.5.8)
    {
#pragma omp task shared(iterator, blocks)
        merge_blocks_recursive(iterator, blocks, blocks_begin, blocks_mid, compare);
#pragma omp task shared(iterator, blocks)
        merge_blocks_recursive(iterator, blocks, blocks_mid, blocks_end, compare);
//#pragma omp taskyield
#pragma omp taskwait
    }
    auto begin = iterator + blocks[blocks_begin];
    auto mid   = iterator + blocks[blocks_mid];
    auto end   = iterator + blocks[blocks_end];
    std::inplace_merge(begin, mid, end, compare);
}
#endif

template <typename RandomAccessIterator, typename Indexable, typename Compare>
void merge_blocks_recursive_seq(RandomAccessIterator& iterator, const Indexable& blocks, size_t blocks_begin,
                                size_t blocks_end, Compare compare) {
    if (blocks_end <= blocks_begin + 1) {
        // recursion done, go back out
        return;
    }
    size_t blocks_mid = (blocks_begin + blocks_end) / 2;
    {
        merge_blocks_recursive_seq(iterator, blocks, blocks_begin, blocks_mid, compare);
        merge_blocks_recursive_seq(iterator, blocks, blocks_mid, blocks_end, compare);
    }
    auto begin = iterator + blocks[blocks_begin];
    auto mid   = iterator + blocks[blocks_mid];
    auto end   = iterator + blocks[blocks_end];

    std::inplace_merge(begin, mid, end, compare);
}

}  // namespace detail

template <typename RandomAccessIterator, typename Compare>
void sort(RandomAccessIterator first, RandomAccessIterator last, Compare compare) {
#if ATLAS_HAVE_OMP_SORTING
    if (atlas_omp_get_max_threads() > 1) {
#pragma omp parallel
#pragma omp single
        detail::merge_sort_recursive(first, 0, std::distance(first, last), compare);
    }
    else {
        std::sort(first, last, compare);
    }
#else
    std::sort(first, last, compare);
#endif
}

template <typename RandomAccessIterator>
void sort(RandomAccessIterator first, RandomAccessIterator last) {
    using value_type = typename std::iterator_traits<RandomAccessIterator>::value_type;
    ::atlas::omp::sort(first, last, std::less<value_type>());
}

template <typename RandomAccessIterator, typename RandomAccessIterator2, typename Compare>
void merge_blocks(RandomAccessIterator first, RandomAccessIterator last, RandomAccessIterator2 blocks_size_first,
                  RandomAccessIterator2 blocks_size_last, Compare compare) {
    using size_type     = typename std::iterator_traits<RandomAccessIterator2>::value_type;
    size_type nb_blocks = std::distance(blocks_size_first, blocks_size_last);
    std::vector<size_type> blocks_displs(nb_blocks + 1);
    blocks_displs[0] = 0;
    for (size_t i = 1; i < blocks_displs.size(); ++i) {
        blocks_displs[i] = blocks_displs[i - 1] + blocks_size_first[i - 1];
    }
#if ATLAS_HAVE_OMP_SORTING
    if (atlas_omp_get_max_threads() > 1) {
#pragma omp parallel
#pragma omp single
        detail::merge_blocks_recursive(first, blocks_displs, 0, nb_blocks, compare);
    }
    else {
        detail::merge_blocks_recursive_seq(first, blocks_displs, 0, nb_blocks, compare);
    }
#else
    detail::merge_blocks_recursive_seq(first, blocks_displs, 0, nb_blocks, compare);
#endif
}

template <typename RandomAccessIterator, typename RandomAccessIterator2>
void merge_blocks(RandomAccessIterator first, RandomAccessIterator last, RandomAccessIterator2 blocks_size_first,
                  RandomAccessIterator2 blocks_size_last) {
    using value_type = typename std::iterator_traits<RandomAccessIterator>::value_type;
    ::atlas::omp::merge_blocks(first, last, blocks_size_first, blocks_size_last, std::less<value_type>());
}

}  // namespace omp
}  // namespace atlas

#undef ATLAS_OMP_TASK_UNTIED_SUPPORTED
#undef ATLAS_HAVE_OMP_SORTING
