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
#include <initializer_list>
#include <vector>
#include "atlas/library/config.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

class ArrayStrides : public std::vector<idx_t> {
private:
    using Base = std::vector<idx_t>;

public:
    ArrayStrides() {}
    ArrayStrides(std::initializer_list<idx_t> list): Base(list) {}
    ArrayStrides(Base&& base): Base(std::forward<Base>(base)) {}
};

namespace detail {

template <typename Int>
inline ArrayStrides make_strides(Int size1) {
    return ArrayStrides{static_cast<idx_t>(size1)};
}
template <typename Int1, typename Int2>
inline ArrayStrides make_strides(Int1 size1, Int2 size2) {
    return ArrayStrides{static_cast<idx_t>(size1), static_cast<idx_t>(size2)};
}
template <typename Int1, typename Int2, typename Int3>
inline ArrayStrides make_strides(Int1 size1, Int2 size2, Int3 size3) {
    return ArrayStrides{static_cast<idx_t>(size1), static_cast<idx_t>(size2), static_cast<idx_t>(size3)};
}
template <typename Int1, typename Int2, typename Int3, typename Int4>
inline ArrayStrides make_strides(Int1 size1, Int2 size2, Int3 size3, Int4 size4) {
    return ArrayStrides{static_cast<idx_t>(size1), static_cast<idx_t>(size2), static_cast<idx_t>(size3),
                        static_cast<idx_t>(size4)};
}
template <typename Int1, typename Int2, typename Int3, typename Int4, typename Int5>
inline ArrayStrides make_strides(Int1 size1, Int2 size2, Int3 size3, Int4 size4, Int5 size5) {
    return ArrayStrides{static_cast<idx_t>(size1), static_cast<idx_t>(size2), static_cast<idx_t>(size3),
                        static_cast<idx_t>(size4), static_cast<idx_t>(size5)};
}

}  // namespace detail

inline ArrayStrides make_strides(std::initializer_list<idx_t> sizes) {
    return ArrayStrides(sizes);
}

template <typename... idx_t>
ArrayStrides make_strides(idx_t... indices) {
    return detail::make_strides(std::forward<idx_t>(indices)...);
}


//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
