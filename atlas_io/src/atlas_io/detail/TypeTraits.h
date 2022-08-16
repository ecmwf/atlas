/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <type_traits>


namespace atlas {
namespace io {

class Metadata;
class Data;
class Stream;

namespace adl_tests {
using std::declval;
using std::is_pod;

template <class, class>
std::false_type can_interprete(...) noexcept(false);

template <class T, class A, class = decltype(interprete(declval<const T&>(), declval<A&>()))>
std::true_type can_interprete(int) noexcept(noexcept(interprete(declval<const T&>(), declval<A&>())));

template <class>
std::false_type can_encode_data(...) noexcept(false);

template <class T, class = decltype(encode_data(declval<const T&>(), declval<Data&>()))>
std::true_type can_encode_data(int) noexcept(noexcept(encode_data(declval<const T&>(), declval<Data&>())));

template <class>
std::false_type can_encode_metadata(...) noexcept(false);

template <class T, class = decltype(encode_metadata(declval<const T&>(), declval<atlas::io::Metadata&>()))>
std::true_type can_encode_metadata(int) noexcept(noexcept(encode_metadata(declval<const T&>(),
                                                                          declval<atlas::io::Metadata&>())));


template <class>
std::false_type can_decode(...) noexcept(false);

template <class T, class = decltype(decode(declval<const Metadata&>(), declval<const Data&>(), declval<T&>()))>
std::true_type can_decode(int) noexcept(noexcept(decode(declval<const Metadata&>(), declval<const Data&>(),
                                                        declval<T&>())));
}  // namespace adl_tests


template <class T, class A>
static constexpr bool is_interpretable() {
    return decltype(adl_tests::can_interprete<typename std::decay<T>::type, A>(0))::value;
}


template <class T>
static constexpr bool can_encode_data() {
    return decltype(adl_tests::can_encode_data<typename std::decay<T>::type>(0))::value;
}

template <class T>
static constexpr bool can_encode_metadata() {
    return decltype(adl_tests::can_encode_metadata<typename std::decay<T>::type>(0))::value;
}


template <typename T>
static constexpr bool is_encodable() {
    return can_encode_data<T>() && can_encode_metadata<T>();
}

template <typename T>
static constexpr bool is_decodable() {
    return decltype(adl_tests::can_decode<typename std::decay<T>::type>(0))::value;
}


template <typename T>
static constexpr bool is_encodable_rvalue() {
    return is_encodable<T>() && std::is_rvalue_reference<T&&>::value;
};

template <bool pred>
using enable_if_t = typename std::enable_if<pred, int>::type;

template <bool pred>
using disable_if_t = typename std::enable_if<!pred, int>::type;

template <typename T>
using enable_if_encodable_t = enable_if_t<is_encodable<T>()>;

template <typename T>
using disable_if_encodable_t = disable_if_t<is_encodable<T>()>;

template <typename T>
using enable_if_decodable_t = enable_if_t<is_decodable<T>()>;

template <typename T>
using disable_if_decodable_t = disable_if_t<is_decodable<T>()>;

template <typename T>
using enable_if_can_encode_metadata_t = enable_if_t<can_encode_metadata<T>()>;

template <typename T>
using disable_if_can_encode_metadata_t = disable_if_t<can_encode_metadata<T>()>;

template <typename T>
using enable_if_can_encode_data_t = enable_if_t<can_encode_data<T>()>;

template <typename T>
using disable_if_can_encode_data_t = disable_if_t<can_encode_data<T>()>;

template <typename T, typename A>
using enable_if_interpretable_t = enable_if_t<is_interpretable<T, A>()>;

template <typename T, typename A>
using disable_if_interpretable_t = disable_if_t<is_interpretable<T, A>()>;

template <typename T>
using enable_if_rvalue_t = enable_if_t<std::is_rvalue_reference<T>::value>;


template <typename T>
using enable_if_move_constructible_encodable_rvalue_t =
    enable_if_t<is_encodable<T>() && std::is_rvalue_reference<T&&>() && std::is_move_constructible<T>()>;

template <typename T>
using enable_if_move_constructible_decodable_rvalue_t =
    enable_if_t<is_decodable<T>() && std::is_rvalue_reference<T&&>() && std::is_move_constructible<T>()>;


}  // namespace io
}  // namespace atlas
