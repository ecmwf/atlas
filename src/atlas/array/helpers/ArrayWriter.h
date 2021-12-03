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

#include "atlas/array.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace helpers {

//------------------------------------------------------------------------------

/// Helper to assign a value to an array or arrayview
// template <typename Value>
struct array_writer;

//------------------------------------------------------------------------------

// Recursive function to apply value to every index
template <typename Value, int Rank, int Dim>
struct array_writer_impl {
    template <typename View, typename... DimIndex>
    static void apply(const View& arr, std::ostream& out, DimIndex... idxs) {
        for (idx_t i = 0; i < arr.shape(Dim); ++i) {
            array_writer_impl<Value, Rank, Dim + 1>::apply(arr, out, idxs..., i);
            if (i < arr.shape(Dim) - 1)
                out << " ";
        }
    }
};

// End of recursion when Dim == Rank
template <typename Value, int Rank>
struct array_writer_impl<Value, Rank, Rank> {
    template <typename View, typename... DimIndex>
    static void apply(const View& arr, std::ostream& out, DimIndex... idxs) {
        out << arr(idxs...);
    }
};

//------------------------------------------------------------------------------

struct array_writer {
    template <typename Value, int Rank>
    static void apply(const ArrayView<Value, Rank>& arr, std::ostream& out) {
        array_writer_impl<Value, Rank, 0u>::apply(arr, out);
        // Note: no need to apply variadic pack (idxs...)
    }

    template <typename Value, int Rank>
    static void apply(const LocalView<Value, Rank>& arr, std::ostream& out) {
        array_writer_impl<Value, Rank, 0u>::apply(arr, out);
        // Note: no need to apply variadic pack (idxs...)
    }
};

//------------------------------------------------------------------------------

}  // namespace helpers
}  // namespace array
}  // namespace atlas
