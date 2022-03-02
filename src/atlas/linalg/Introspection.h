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
#include "eckit/linalg/types.h"

#include "atlas/library/config.h"

namespace atlas {
namespace linalg {
namespace introspection {


namespace detail {
// void_t is a C++17 feature
template <typename... Ts>
struct make_void {
    typedef void type;
};
template <typename... Ts>
using void_t = typename make_void<Ts...>::type;
}  // namespace detail

template <typename, typename = void>
struct has_contiguous : std::false_type {};

template <typename T>
struct has_contiguous<T, detail::void_t<decltype(&T::contiguous)>> : std::true_type {};

template <typename, typename = void>
struct has_rank : std::false_type {};

template <typename T>
struct has_rank<T, detail::void_t<decltype(&T::rank)>> : std::true_type {};

template <typename, typename = void>
struct has_RANK : std::false_type {};

template <typename T>
struct has_RANK<T, detail::void_t<decltype(&T::RANK)>> : std::true_type {};

template <typename, typename = void>
struct has_shape : std::false_type {};

template <typename T>
struct has_shape<T, detail::void_t<decltype(std::declval<T>().shape(0))>> : std::true_type {};

namespace detail {
template <typename View>
struct BaseType {
    template <typename Base>
    static constexpr int is_derived_from() {
        return std::is_base_of<Base, View>::value;
    }

    using type = typename std::conditional<
        is_derived_from<eckit::linalg::Matrix>(), eckit::linalg::Matrix,
        typename std::conditional<is_derived_from<eckit::linalg::Vector>(), eckit::linalg::Vector,
                                  typename std::remove_const<View>::type /*default*/
                                  >::type>::type;
};
}  // namespace detail

template <typename View>
using base_type = typename detail::BaseType<View>::type;

namespace detail {
template <typename T, typename>
struct ValueType {
    using type = typename T::value_type;
};

template <typename T>
struct ValueType<T, eckit::linalg::Matrix> {
    using type = eckit::linalg::Scalar;
};

template <typename T>
struct ValueType<T, eckit::linalg::Vector> {
    using type = eckit::linalg::Scalar;
};

template <typename View, typename>
struct Introspect {
    template <typename T, typename std::enable_if<has_contiguous<T>::value, T>::type* = nullptr>
    static bool contiguous(const T& view) {
        return view.contiguous();
    }
    template <typename T, typename std::enable_if<not has_contiguous<T>::value, T>::type* = nullptr>
    static bool contiguous(const T&) {
        return true;
    }

    template <typename T, typename std::enable_if<has_rank<T>::value, T>::type* = nullptr>
    static constexpr int rank() {
        static_assert(has_RANK<T>::value, "workaround for Cray 8.7: we need to use T::RANK instead of T::rank()");
        return T::RANK;
    }
    template <typename T, typename std::enable_if<not has_rank<T>::value, T>::type* = nullptr>
    static constexpr int rank() {
        return 1;
    }

    static constexpr int rank() { return rank<View>(); }

    template <typename T, typename std::enable_if<has_shape<T>::value, T>::type* = nullptr>
    static idx_t outer_dimension(const T& view) {
        return view.shape(0);
    }
    template <typename T, typename std::enable_if<not has_shape<T>::value && (rank<T>() == 1), T>::type* = nullptr>
    static idx_t outer_dimension(const T& view) {
        return view.size();
    }

    template <typename T, typename std::enable_if<has_shape<T>::value, T>::type* = nullptr>
    static idx_t inner_dimension(const T& view) {
        return view.shape(rank<T>() - 1);
    }
    template <typename T, typename std::enable_if<not has_shape<T>::value && (rank<T>() == 1), T>::type* = nullptr>
    static idx_t inner_dimension(const T& view) {
        return view.size();
    }

    template <int Dim, typename T, typename std::enable_if<not has_shape<T>::value, T>::type* = nullptr>
    static idx_t shape(const T& view) {
        //static_assert( Dim == 0, "Dim must be 0 if shape() is not there" );
        return view.size();
    }

    template <int Dim, typename T, typename std::enable_if<has_shape<T>::value, T>::type* = nullptr>
    static idx_t shape(const T& view) {
        return view.shape(Dim);
    }

    static constexpr bool layout_left() { return true; }

    using value_type = typename ValueType<View, base_type<View>>::type;
};

template <typename View>
struct Introspect<View, eckit::linalg::Matrix> {
    using value_type = eckit::linalg::Scalar;
    static constexpr int rank() { return 2; }
    static bool contiguous(const eckit::linalg::Matrix&) { return true; }
    static idx_t inner_dimension(const eckit::linalg::Matrix& m) { return static_cast<idx_t>(m.cols()); }
    static idx_t outer_dimension(const eckit::linalg::Matrix& m) { return static_cast<idx_t>(m.rows()); }
    template <int Dim>
    static idx_t shape(const eckit::linalg::Matrix& m) {
        static_assert(Dim <= 1, "Dim must be 0 or 1");
        return (Dim == 0) ? outer_dimension(m) : inner_dimension(m);
    }
    static constexpr bool layout_left() { return false; }
};

template <typename T>
using Introspection = Introspect<T, base_type<T>>;

}  // namespace detail

template <typename T>
using value_type = typename detail::Introspection<T>::value_type;

template <typename View>
using base_type = typename detail::BaseType<View>::type;

template <typename T>
static bool contiguous(const T& view) {
    return detail::Introspection<T>::contiguous(view);
}
template <typename T>
static constexpr int rank() {
    return detail::Introspection<T>::rank();
}
template <typename T>
static idx_t outer_dimension(const T& view) {
    return detail::Introspection<T>::outer_dimension(view);
}
template <typename T>
static idx_t inner_dimension(const T& view) {
    return detail::Introspection<T>::inner_dimension(view);
}

template <int Dim, typename T>
static idx_t shape(const T& view) {
    return detail::Introspection<T>::template shape<Dim>(view);
}

template <typename T>
static constexpr bool layout_left(const T&) {
    return detail::Introspection<T>::layout_left();
}
template <typename T>
static constexpr bool layout_right(const T&) {
    return not detail::Introspection<T>::layout_left();
}

}  // namespace introspection
}  // namespace linalg
}  // namespace atlas
