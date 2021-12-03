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

#include <stddef.h>
#include <utility>
#include <vector>

#include "atlas/library/config.h"

//------------------------------------------------------------------------------------------------------


namespace atlas {
namespace array {

//using ArrayIdx = std::vector<idx_t>;
typedef std::vector<idx_t> ArrayIdx;

namespace detail {

inline ArrayIdx make_idx(idx_t size1) {
    return std::vector<idx_t>(1, size1);
}
inline ArrayIdx make_idx(idx_t size1, idx_t size2) {
    std::vector<idx_t> v(2);
    v[0] = size1;
    v[1] = size2;
    return v;
}
inline ArrayIdx make_idx(idx_t size1, idx_t size2, idx_t size3) {
    std::vector<idx_t> v(3);
    v[0] = size1;
    v[1] = size2;
    v[2] = size3;
    return v;
}
inline ArrayIdx make_idx(idx_t size1, idx_t size2, idx_t size3, idx_t size4) {
    std::vector<idx_t> v(4);
    v[0] = size1;
    v[1] = size2;
    v[2] = size3;
    v[3] = size4;
    return v;
}

}  // namespace detail

template <typename... idx_t>
ArrayIdx make_idx(idx_t... indices) {
    return detail::make_idx(std::forward<idx_t>(indices)...);
}


//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
