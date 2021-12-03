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
#include <iterator>

#include "atlas/parallel/omp/omp.h"

namespace atlas {
namespace omp {

template <typename BidirIt, typename T>
void fill(BidirIt begin, BidirIt end, const T& value) {
    if (atlas_omp_get_max_threads() > 1) {
        auto size = std::distance(begin, end);
        atlas_omp_parallel {
            auto nthreads  = atlas_omp_get_num_threads();
            auto tid       = atlas_omp_get_thread_num();
            auto chunksize = size / nthreads;
            auto v_begin   = begin + chunksize * tid;
            auto v_end     = (tid == nthreads - 1) ? end : v_begin + chunksize;
            std::fill(v_begin, v_end, value);
        }
    }
    else {
        std::fill(begin, end, value);
    }
}

}  // namespace omp
}  // namespace atlas
