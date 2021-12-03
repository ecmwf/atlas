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

#include <cstddef>
#include <initializer_list>
#include <vector>

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace array {

class ArrayLayout : public std::vector<idx_t> {
private:
    using Base = std::vector<idx_t>;

public:
    ArrayLayout() {}
    ArrayLayout(std::initializer_list<idx_t> list): Base(list) {}
    ArrayLayout(Base&& base): Base(std::forward<Base>(base)) {}
};

namespace detail {

inline ArrayLayout make_layout(idx_t size1) {
    return ArrayLayout{size1};
}
inline ArrayLayout make_layout(idx_t size1, idx_t size2) {
    return ArrayLayout{size1, size2};
}
inline ArrayLayout make_layout(idx_t size1, idx_t size2, idx_t size3) {
    return ArrayLayout{size1, size2, size3};
}
inline ArrayLayout make_layout(idx_t size1, idx_t size2, idx_t size3, idx_t size4) {
    return ArrayLayout{size1, size2, size3, size4};
}
inline ArrayLayout make_layout(idx_t size1, idx_t size2, idx_t size3, idx_t size4, idx_t size5) {
    return ArrayLayout{size1, size2, size3, size4, size5};
}

}  // namespace detail

template <typename... idx_t>
ArrayLayout make_layout(idx_t... indices) {
    return detail::make_layout(std::forward<idx_t>(indices)...);
}


//------------------------------------------------------------------------------------------------------

}  // namespace array
}  // namespace atlas
