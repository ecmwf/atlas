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

#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/LocalView.h"
#include "atlas/array_fwd.h"

namespace atlas {
namespace array {

extern template IndexView<idx_t, 1> make_indexview<idx_t, 1>(Array&);
extern template IndexView<idx_t, 2> make_indexview<idx_t, 2>(Array&);
extern template IndexView<const idx_t, 1> make_indexview<const idx_t, 1>(Array&);
extern template IndexView<const idx_t, 2> make_indexview<const idx_t, 2>(Array&);
extern template IndexView<const idx_t, 1> make_indexview<idx_t, 1>(const Array&);
extern template IndexView<const idx_t, 2> make_indexview<idx_t, 2>(const Array&);

#define EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK(TYPE, RANK)                                                           \
    extern template ArrayView<TYPE, RANK> make_view<TYPE, RANK>(Array&);                                              \
    extern template ArrayView<const TYPE, RANK> make_view<const TYPE, RANK>(Array&);                                  \
    extern template ArrayView<const TYPE, RANK> make_view<TYPE, RANK>(const Array&);                                  \
    extern template ArrayView<const TYPE, RANK> make_view<const TYPE, RANK>(const Array&);                            \
                                                                                                                      \
    extern template LocalView<TYPE, RANK> make_view<TYPE, RANK, nullptr>(TYPE data[], const ArrayShape&);             \
    extern template LocalView<const TYPE, RANK> make_view<const TYPE, RANK, nullptr>(TYPE data[], const ArrayShape&); \
    extern template LocalView<const TYPE, RANK> make_view<TYPE, RANK, nullptr>(const TYPE data[], const ArrayShape&); \
    extern template LocalView<const TYPE, RANK> make_view<const TYPE, RANK, nullptr>(const TYPE data[],               \
                                                                                     const ArrayShape&);              \
                                                                                                                      \
    extern template LocalView<TYPE, RANK> make_view<TYPE, RANK, nullptr>(TYPE data[], size_t);                        \
    extern template LocalView<const TYPE, RANK> make_view<const TYPE, RANK, nullptr>(TYPE data[], size_t);            \
    extern template LocalView<const TYPE, RANK> make_view<TYPE, RANK, nullptr>(const TYPE data[], size_t);            \
    extern template LocalView<const TYPE, RANK> make_view<const TYPE, RANK, nullptr>(const TYPE data[], size_t);


#define EXPLICIT_TEMPLATE_DECLARATION(RANK)              \
    EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK(int, RANK)   \
    EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK(long, RANK)  \
    EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK(float, RANK) \
    EXPLICIT_TEMPLATE_DECLARATION_TYPE_RANK(double, RANK)


// For each NDims in [1..9]
EXPLICIT_TEMPLATE_DECLARATION(1)
EXPLICIT_TEMPLATE_DECLARATION(2)
EXPLICIT_TEMPLATE_DECLARATION(3)
EXPLICIT_TEMPLATE_DECLARATION(4)
EXPLICIT_TEMPLATE_DECLARATION(5)
EXPLICIT_TEMPLATE_DECLARATION(6)
EXPLICIT_TEMPLATE_DECLARATION(7)
EXPLICIT_TEMPLATE_DECLARATION(8)
EXPLICIT_TEMPLATE_DECLARATION(9)

#undef EXPLICIT_TEMPLATE_DECLARATION

}  // namespace array
}  // namespace atlas
