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

#include "eckit/linalg/Matrix.h"
#include "eckit/linalg/Vector.h"

#include "atlas/array/LocalView.h"
#include "atlas/linalg/Introspection.h"

namespace atlas {
namespace linalg {

template <typename Value, int Rank>
using View = array::LocalView<Value, Rank>;

namespace view {

template <typename View, typename>
struct Convert {
    using value_type = typename std::conditional<std::is_const<View>::value, const introspection::value_type<View>,
                                                 introspection::value_type<View>>::type;
    using type       = array::LocalView<value_type, introspection::rank<View>()>;
    static type apply(View& view) { return type(view.data(), view.shape()); }
};

template <typename View>
struct Convert<View, eckit::linalg::Vector> {
    using value_type =
        typename std::conditional<std::is_const<View>::value, const eckit::linalg::Scalar, eckit::linalg::Scalar>::type;
    using type = array::LocalView<value_type, 1>;
    static type apply(View& v) { return type(v.data(), array::make_shape(v.size())); }
};

template <typename View>
struct Convert<View, eckit::linalg::Matrix> {
    using value_type =
        typename std::conditional<std::is_const<View>::value, const eckit::linalg::Scalar, eckit::linalg::Scalar>::type;
    using type = array::LocalView<value_type, 2>;
    static type apply(View& m) { return type(m.data(), array::make_shape(m.cols(), m.rows())); }
};

template <typename View>
struct ConvertView {
    using type = typename Convert<View, introspection::base_type<View>>::type;
    static type apply(View& v) { return Convert<View, introspection::base_type<View>>::apply(v); }
};

}  // namespace view

template <typename View>
using view_type = typename view::ConvertView<View>::type;

template <typename T>
view_type<T> make_view(T& view) {
    return view::ConvertView<T>::apply(view);
}

}  // namespace linalg
}  // namespace atlas
