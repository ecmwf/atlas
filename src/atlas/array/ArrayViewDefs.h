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

#include <utility>

namespace atlas {
namespace array {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <int cDim>
struct Dim {
    static constexpr int cdim = cDim;
};

struct LastDim {};
struct FirstDim {};

template <typename T>
struct is_dim_policy : std::false_type {};
template <int cDim>
struct is_dim_policy<Dim<cDim>> : std::true_type {};

template <>
struct is_dim_policy<LastDim> : std::true_type {};
template <>
struct is_dim_policy<FirstDim> : std::true_type {};
#endif
}  // namespace array
}  // namespace atlas
