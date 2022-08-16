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

#ifndef ATLAS_IO_STATIC_ASSERT
#if defined(__clang_analyzer__)
#define ATLAS_IO_STATIC_ASSERT 0
#elif defined(__clang__) && (__clang_major__ < 4)
#define ATLAS_IO_STATIC_ASSERT 0
#else
#define ATLAS_IO_STATIC_ASSERT 1
#endif
#endif

#if ATLAS_IO_STATIC_ASSERT

#include <cstdlib>
#include "atlas_io/detail/TypeTraits.h"

namespace atlas {
namespace io {

template <typename T, disable_if_can_encode_data_t<T> = 0>
void encode_data(const T& in, atlas::io::Data& out) {
    static_assert(
        can_encode_data<T>(),
        "\n\n"
        "\n     Static assertion failed"
        "\n     -----------------------"
        "\n"
        "\n     Values of template type T cannot be encoded into Data because following function is not defined:"
        "\n"
        "\n         void encode_data(const T& value, atlas::io::Data& out);"
        "\n"
        "\n     Note that argument-dependent-lookup rules apply."
        "\n     --> The function must be declared in the namespace of type T."
        "\n\n");
}


template <typename T, disable_if_can_encode_metadata_t<T> = 0>
size_t encode_metadata(const T&, atlas::io::Metadata&) {
    static_assert(
        can_encode_metadata<T>(),
        "\n\n"
        "\n     Static assertion failed"
        "\n     -----------------------"
        "\n"
        "\n     Values of template type T cannot be incoded into Metadata because following function is not defined:"
        "\n"
        "\n         size_t encode_metadata(const T& value, atlas::io::Metadata& metadata);"
        "\n"
        "\n     Note that argument-dependent-lookup rules apply."
        "\n     --> The function must be declared in the namespace of type T"
        "\n\n");
    return 0;
}


template <typename T, disable_if_decodable_t<T> = 0>
void decode(const atlas::io::Metadata&, const atlas::io::Data&, T&) {
    static_assert(is_decodable<T>(),
                  "\n\n"
                  "\n     Static assertion failed"
                  "\n     -----------------------"
                  "\n"
                  "\n     Values of template type T cannot be decoded because following function is not defined:"
                  "\n"
                  "\n         void decode(const atlas::io::Metadata&, const atlas::io::Data&, T& out);"
                  "\n"
                  "\n     Note that argument-dependent-lookup rules apply."
                  "\n     The function must be declared in the namespace of type T"
                  "\n\n");
}

}  // namespace io
}  // namespace atlas

#endif
