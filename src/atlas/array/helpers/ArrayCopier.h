/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/array.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace helpers {

//------------------------------------------------------------------------------

/// Helper to copy from source array to target array.
template <typename Value, unsigned int Rank>
struct array_copier;

//------------------------------------------------------------------------------

// Recursive function to copy array elements across all dimensions.
template <typename Value, unsigned int Rank, unsigned int Dim>
struct array_copier_impl {
    template <typename SourceView, typename TargetView, typename... DimIndex>
    static void apply(const SourceView& sourceArr, TargetView& targetArr, const std::array<idx_t, Rank>& shape,
                      DimIndex... idxs) {
        for (idx_t i = 0; i < shape[Dim]; ++i) {
            array_copier_impl<Value, Rank, Dim + 1>::apply(sourceArr, targetArr, shape, idxs..., i);
        }
    }
};

// End of recursion when Dim == Rank
template <typename Value, unsigned int Rank>
struct array_copier_impl<Value, Rank, Rank> {
    template <typename SourceView, typename TargetView, typename... DimIndex>
    static void apply(const SourceView& sourceArr, TargetView& targetArr, const std::array<idx_t, Rank>& shape,
                      DimIndex... idxs) {
        targetArr(idxs...) = sourceArr(idxs...);
    }
};

//------------------------------------------------------------------------------

template <typename Value, unsigned int Rank>
struct array_copier {
    // Copy from source array to target array.
    static void apply(const ArrayView<Value, Rank>& sourceArr, ArrayView<Value, Rank>& targetArr) {
        array_copier_impl<Value, Rank, 0u>::apply(sourceArr, targetArr, shape(sourceArr, targetArr));
    }
    // Copy from const source array to target array.
    static void apply(const ArrayView<const Value, Rank>& sourceArr, ArrayView<Value, Rank>& targetArr) {
        array_copier_impl<Value, Rank, 0u>::apply(sourceArr, targetArr, shape(sourceArr, targetArr));
    }

private:
    // Make an array shape that is compatible with sourceArr and targetArr.
    // This is useful for arrays from fields with different sized halos.
    template <typename SourceView, typename TargetView>
    static std::array<idx_t, Rank> shape(const SourceView& sourceArr, const TargetView& targetArr) {
        auto arrShape = std::array<idx_t, Rank>{};
        for (unsigned int Dim = 0; Dim < Rank; ++Dim) {
            arrShape[Dim] = std::min(sourceArr.shape(Dim), targetArr.shape(Dim));
        }
        return arrShape;
    }
};

//------------------------------------------------------------------------------

}  // namespace helpers
}  // namespace array
}  // namespace atlas
