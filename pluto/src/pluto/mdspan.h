/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <type_traits>

#include "pluto/pluto_config.h"

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wpre-c++2b-compat"
#endif

#if PLUTO_HAVE_MDSPAN
#include <mdspan>

#pragma push_macro("STD_MDSPAN_NAMESPACE")
#undef STD_MDSPAN_NAMESPACE
#define STD_MDSPAN_NAMESPACE std
namespace pluto {
using ::STD_MDSPAN_NAMESPACE::dynamic_extent;
using ::STD_MDSPAN_NAMESPACE::layout_left;
using ::STD_MDSPAN_NAMESPACE::layout_right;
using ::STD_MDSPAN_NAMESPACE::layout_stride;
using ::STD_MDSPAN_NAMESPACE::default_accessor;
using ::STD_MDSPAN_NAMESPACE::extents;
using ::STD_MDSPAN_NAMESPACE::dextents;

#if PLUTO_MDSPAN_USE_PAREN_OPERATOR
#include "pluto/detail/mdspan_paren_operator.h"
#else
using ::STD_MDSPAN_NAMESPACE::mdspan;
#endif
} // namespace pluto
#pragma pop_macro("STD_MDSPAN_NAMESPACE")
#define PLUTO_MDSPAN_USE_BRACKET_OPERATOR 1
#define PLUTO_MDSPAN_HOST_DEVICE

#else

#pragma push_macro("MDSPAN_IMPL_STANDARD_NAMESPACE")
#undef MDSPAN_IMPL_STANDARD_NAMESPACE
#define MDSPAN_IMPL_STANDARD_NAMESPACE pluto
#include "pluto/detail/mdspan/mdspan.hpp"
#pragma pop_macro("MDSPAN_IMPL_STANDARD_NAMESPACE")
#define PLUTO_MDSPAN_USE_BRACKET_OPERATOR MDSPAN_USE_BRACKET_OPERATOR
#define PLUTO_MDSPAN_HOST_DEVICE MDSPAN_IMPL_HOST_DEVICE
#endif

namespace pluto {
    // A C++26 addition:
    template< std::size_t Rank, class IndexType = std::size_t >
    using dims = dextents<IndexType, Rank>;

    // A C++26 addition:
    template<std::size_t ByteAlignment, class ElementType>
    bool is_sufficiently_aligned(ElementType* ptr) {
        static_assert(ByteAlignment != 0 && (ByteAlignment & (ByteAlignment - 1)) == 0,
                      "ByteAlignment must be a power of two.");
        return reinterpret_cast<std::uintptr_t>(ptr) % ByteAlignment == 0;
    }

    // A C++26 addition:
    template<class ElementType, std::size_t ByteAlignment>
    struct aligned_accessor {
        static_assert(ByteAlignment != 0 && (ByteAlignment & (ByteAlignment - 1)) == 0,
                      "ByteAlignment must be a power of two.");
        static_assert(ByteAlignment >= alignof(ElementType), "Insufficient byte alignment for ElementType.");

        using offset_policy    = default_accessor<ElementType>;
        using element_type     = ElementType;
        using reference        = ElementType&;
        using data_handle_type = ElementType*;

        static constexpr std::size_t byte_alignment = ByteAlignment;

        constexpr aligned_accessor() noexcept = default;

        template<
            class OtherElementType,
            std::size_t OtherByteAlignment,
            typename = std::enable_if_t<
                std::is_convertible<OtherElementType(*)[], element_type(*)[]>::value &&
                (OtherByteAlignment >= byte_alignment)> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr aligned_accessor(aligned_accessor<OtherElementType, OtherByteAlignment>) noexcept {}

        template<
            class OtherElementType,
            typename = std::enable_if_t<std::is_convertible<OtherElementType(*)[], element_type(*)[]>::value> >
        PLUTO_MDSPAN_HOST_DEVICE
        explicit constexpr aligned_accessor(default_accessor<OtherElementType>) noexcept {}

        template<
            class OtherElementType,
            typename = std::enable_if_t<std::is_convertible<element_type(*)[], OtherElementType(*)[]>::value> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr operator default_accessor<OtherElementType>() const noexcept {
            return {};
        }

        PLUTO_MDSPAN_HOST_DEVICE
        constexpr reference access(data_handle_type ptr, std::size_t index) const noexcept {
            return assume_aligned(ptr)[index];
        }

        PLUTO_MDSPAN_HOST_DEVICE
        constexpr typename offset_policy::data_handle_type offset(data_handle_type ptr, std::size_t index) const noexcept {
            return assume_aligned(ptr) + index;
        }

    private:
        PLUTO_MDSPAN_HOST_DEVICE
        static inline constexpr auto assume_aligned(data_handle_type ptr) noexcept {
        #if defined(__cpp_lib_assume_aligned)
            return std::assume_aligned<byte_alignment>(ptr);
        #elif defined(__GNUC__) || defined(__clang__)
            return static_cast<data_handle_type>(__builtin_assume_aligned(ptr, byte_alignment));
        #else
            return ptr;
        #endif
        }
    };

// Platform macro for cross-compiler restrict support
#if defined(__GNUC__) || defined(__clang__)
    #define PLUTO_MDSPAN_RESTRICT __restrict__
#elif defined(__INTEL_COMPILER)
    #define PLUTO_MDSPAN_RESTRICT __restrict
#else
    #define PLUTO_MDSPAN_RESTRICT
#endif

    // This is **not** part of the C++ standard, but is a common extension to mdspan implementations.
    // The restrict keyword is not standard C++, but is supported by many compilers.
    //It indicates that the pointer is the only reference to the object it points to, which can enable certain optimizations.
    template<class ElementType>
    struct restrict_accessor {
        using element_type     = ElementType;
        using reference        = ElementType&;
        using data_handle_type = ElementType* PLUTO_MDSPAN_RESTRICT;
        using offset_policy    = restrict_accessor<ElementType>;

        constexpr restrict_accessor() noexcept = default;

        template<class OtherElementType, typename = std::enable_if_t<std::is_convertible<OtherElementType(*)[], ElementType(*)[]>::value> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr restrict_accessor(restrict_accessor<OtherElementType>) noexcept {}

        template<class OtherElementType, typename = std::enable_if_t<std::is_convertible<OtherElementType(*)[], ElementType(*)[]>::value> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr restrict_accessor(default_accessor<OtherElementType>) noexcept {}

        template<class OtherElementType, typename = std::enable_if_t<std::is_convertible<ElementType(*)[], OtherElementType(*)[]>::value> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr operator default_accessor<OtherElementType>() const noexcept {
            return {};
        }

        PLUTO_MDSPAN_HOST_DEVICE
        constexpr reference access(data_handle_type ptr, std::size_t index) const noexcept {
            return ptr[index];
        }

        PLUTO_MDSPAN_HOST_DEVICE
        constexpr typename offset_policy::data_handle_type offset(data_handle_type ptr, std::size_t index) const noexcept {
            return ptr + index;
        }
    };

// A combination of restrict_accessor and aligned_accessor
    template<class ElementType, std::size_t ByteAlignment>
    struct restrict_aligned_accessor {
        static_assert(ByteAlignment != 0 && (ByteAlignment & (ByteAlignment - 1)) == 0,
                      "ByteAlignment must be a power of two.");
        static_assert(ByteAlignment >= alignof(ElementType), "Insufficient byte alignment for ElementType.");

        using element_type     = ElementType;
        using reference        = ElementType&;
        using data_handle_type = ElementType* PLUTO_MDSPAN_RESTRICT;
        using offset_policy    = restrict_accessor<ElementType>;

        static constexpr std::size_t byte_alignment = ByteAlignment;

        constexpr restrict_aligned_accessor() noexcept = default;

        template<
            class OtherElementType,
            std::size_t OtherByteAlignment,
            typename = std::enable_if_t<
                std::is_convertible<OtherElementType(*)[], ElementType(*)[]>::value &&
                (OtherByteAlignment >= byte_alignment)> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr restrict_aligned_accessor(restrict_aligned_accessor<OtherElementType, OtherByteAlignment>) noexcept {}

        template<
            class OtherElementType,
            typename = std::enable_if_t<std::is_convertible<OtherElementType(*)[], ElementType(*)[]>::value> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr restrict_aligned_accessor(restrict_accessor<OtherElementType>) noexcept {}

        template<
            class OtherElementType,
            std::size_t OtherByteAlignment,
            typename = std::enable_if_t<
                std::is_convertible<OtherElementType(*)[], ElementType(*)[]>::value &&
                (OtherByteAlignment >= byte_alignment)> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr restrict_aligned_accessor(aligned_accessor<OtherElementType, OtherByteAlignment>) noexcept {}

        template<
            class OtherElementType,
            typename = std::enable_if_t<std::is_convertible<OtherElementType(*)[], ElementType(*)[]>::value> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr restrict_aligned_accessor(default_accessor<OtherElementType>) noexcept {}

        template<
            class OtherElementType,
            typename = std::enable_if_t<std::is_convertible<element_type(*)[], OtherElementType(*)[]>::value> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr operator default_accessor<OtherElementType>() const noexcept {
            return {};
        }

        template<
            class OtherElementType,
            typename = std::enable_if_t<std::is_convertible<element_type(*)[], OtherElementType(*)[]>::value> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr operator restrict_accessor<OtherElementType>() const noexcept {
            return {};
        }

        template<
            class OtherElementType,
            typename = std::enable_if_t<std::is_convertible<element_type(*)[], OtherElementType(*)[]>::value> >
        PLUTO_MDSPAN_HOST_DEVICE
        constexpr operator aligned_accessor<OtherElementType, byte_alignment>() const noexcept {
            return {};
        }

        PLUTO_MDSPAN_HOST_DEVICE
        constexpr reference access(data_handle_type ptr, std::size_t index) const noexcept {
            return assume_aligned(ptr)[index];
        }

        PLUTO_MDSPAN_HOST_DEVICE
        constexpr typename offset_policy::data_handle_type offset(data_handle_type ptr, std::size_t index) const noexcept {
            return assume_aligned(ptr) + index;
        }

    private:
        PLUTO_MDSPAN_HOST_DEVICE
        static inline constexpr data_handle_type assume_aligned(data_handle_type ptr) noexcept {
        #if defined(__cpp_lib_assume_aligned)
            return static_cast<data_handle_type>(std::assume_aligned<byte_alignment>(ptr));
        #elif defined(__GNUC__) || defined(__clang__)
            return static_cast<data_handle_type>(__builtin_assume_aligned(ptr, byte_alignment));
        #else
            return ptr;
        #endif
        }
    };
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#if defined(__INTEL_COMPILER)
#pragma warning pop
#endif

// ------------------------------------------------------------------------------------------------
