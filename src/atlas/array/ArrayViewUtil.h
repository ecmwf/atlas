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

namespace atlas {
namespace array {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <int cnt, int DimSkip, typename DATA_TYPE, int RANK>
constexpr typename std::enable_if<(cnt == RANK), idx_t>::type get_var_size_impl(
    array::ArrayView<DATA_TYPE, RANK>& field) {
    return 1;
}

template <int cnt, int DimSkip, typename DATA_TYPE, int RANK>
constexpr typename std::enable_if<(cnt != RANK), idx_t>::type get_var_size_impl(
    array::ArrayView<DATA_TYPE, RANK>& field) {
    return (cnt == DimSkip) ? get_var_size_impl<cnt + 1, DimSkip>(field)
                            : get_var_size_impl<cnt + 1, DimSkip>(field) * field.template shape<cnt>();
}

template <int DimSkip, typename DATA_TYPE, int RANK>
constexpr idx_t get_var_size(array::ArrayView<DATA_TYPE, RANK>& field) {
    return get_var_size_impl<0, DimSkip>(field);
}

template <typename DimPolicy>
struct get_dim {
    static constexpr int value = -1;
};

template <int cDim>
struct get_dim<Dim<cDim>> {
    static constexpr int value = cDim;
};

template <typename DimPolicy, typename DATA_TYPE, int RANK>
constexpr unsigned int get_parallel_dim(array::ArrayView<DATA_TYPE, RANK>& field) {
    static_assert(is_dim_policy<DimPolicy>::value, "DimPolicy required");
    return std::is_same<DimPolicy, FirstDim>::value
               ? 0
               : (std::is_same<DimPolicy, LastDim>::value ? RANK - 1 : (get_dim<DimPolicy>::value));
}
#endif

}  // namespace array
}  // namespace atlas
