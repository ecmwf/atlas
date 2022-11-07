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

#include <type_traits>

#include "atlas/array.h"
#include "atlas/runtime/Exception.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace helpers {

//------------------------------------------------------------------------------

/// Helper to assign a value to an array or arrayview
template <typename Value, unsigned int Rank>
struct array_assigner;

//------------------------------------------------------------------------------

// Recursive function to apply value to every index
template <typename Value, unsigned int Rank, unsigned int Dim>
struct array_assigner_impl {
    template <typename View, typename... DimIndex>
    static void apply(View& arr, Value value, DimIndex... idxs) {
        for (idx_t i = 0; i < arr.shape(Dim); ++i) {
            array_assigner_impl<Value, Rank, Dim + 1>::apply(arr, value, idxs..., i);
        }
    }

    template <typename View, typename Iterator, typename... DimIndex>
    static void apply(View& arr, Iterator& it, DimIndex... idxs) {
        for (idx_t i = 0; i < arr.shape(Dim); ++i) {
            array_assigner_impl<Value, Rank, Dim + 1>::apply(arr, it, idxs..., i);
        }
    }
};

// End of recursion when Dim == Rank
template <typename Value, unsigned int Rank>
struct array_assigner_impl<Value, Rank, Rank> {
    template <typename View, typename... DimIndex>
    static void apply(View& arr, Value value, DimIndex... idxs) {
        arr(idxs...) = value;
    }

    template <typename View, typename Iterator, typename... DimIndex>
    static void apply(View& arr, Iterator& it, DimIndex... idxs) {
        arr(idxs...) = *it;
        ++it;
    }
};

//------------------------------------------------------------------------------

template <typename Value, unsigned int Rank>
struct array_assigner {
    template <typename T>
    static void apply(Array& arr, Value value) {
        return apply(make_host_view<T, Rank>(arr), value);
    }

    static void apply(ArrayView<Value, Rank>& arr, Value value) {
        array_assigner_impl<Value, Rank, 0u>::apply(arr, value);
        // Note: no need to apply variadic pack (idxs...)
    }

    template <typename Iterable>
    static void apply(ArrayView<Value, Rank>& arr, const Iterable& iterable) {
        typename Iterable::const_iterator it = iterable.begin();
        array_assigner_impl<Value, Rank, 0u>::apply(arr, it);
        ATLAS_ASSERT(it = iterable.end());
    }

    template <typename Iterable>
    static void apply(IndexView<Value, Rank>& arr, const Iterable& iterable) {
        typename Iterable::const_iterator it = iterable.begin();
        array_assigner_impl<Value, Rank, 0u>::apply(arr, it);
        ATLAS_ASSERT(it = iterable.end());
    }


    static void apply(LocalView<Value, Rank>& arr, Value value) {
        array_assigner_impl<Value, Rank, 0u>::apply(arr, value);
        // Note: no need to apply variadic pack (idxs...)
    }
};

//------------------------------------------------------------------------------

}  // namespace helpers
}  // namespace array
}  // namespace atlas
