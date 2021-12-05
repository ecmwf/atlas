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

#include <iostream>
#include "atlas/runtime/Log.h"

namespace atlas {
namespace omp {

template <class InputIterator, class OutputIterator>
OutputIterator copy(InputIterator first, InputIterator last, OutputIterator result) {
    if (atlas_omp_get_max_threads() > 1) {
        auto size = std::distance(first, last);
        atlas_omp_parallel {
            auto nthreads     = atlas_omp_get_num_threads();
            auto tid          = atlas_omp_get_thread_num();
            auto chunksize    = size / nthreads;
            auto v_begin      = first + chunksize * tid;
            auto v_end        = (tid == nthreads - 1) ? last : v_begin + chunksize;
            auto result_begin = result + chunksize * tid;
            std::copy(v_begin, v_end, result_begin);
        }
        return result + size;
    }
    else {
        return std::copy(first, last, result);
    }
}

}  // namespace omp
}  // namespace atlas
