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

#include "pluto/mdspan.h"

namespace atlas {
using ::pluto::dynamic_extent;
using ::pluto::layout_left;
using ::pluto::layout_right;
using ::pluto::layout_stride;
using ::pluto::default_accessor;
using ::pluto::extents;
using ::pluto::dextents;
using ::pluto::dims;
using ::pluto::mdspan;
} // namespace atlas

// ------------------------------------------------------------------------------------------------

#include <type_traits>

namespace atlas {

template <typename T, size_t Base>
class index_reference {
public:
    using value_type = T;
    constexpr index_reference(value_type& idx): idx_(idx) {
    }
    constexpr void set(const value_type& value) {
        idx_ = value + base_;
    }
    constexpr value_type get() const {
        return idx_ - base_;
    }
    constexpr void operator =(const value_type& value) {
        set(value);
    }
    constexpr index_reference& operator=(const index_reference& other) noexcept {
        set(other.get());
        return *this;
    }
    constexpr index_reference& operator--() noexcept {
        --idx_;
        return *this;
    }
    constexpr index_reference& operator++() noexcept {
        ++idx_;
        return *this;
    }
    constexpr index_reference& operator+=(value_type v) noexcept {
        idx_ += v;
        return *this;
    }
    constexpr index_reference& operator-=(value_type v) noexcept {
        idx_ -= v;
        return *this;
    }
    constexpr operator value_type() const noexcept{
        return get();
    }

private:
    value_type& idx_;
    static constexpr value_type base_ = Base;
};

template<class ElementType, size_t Base>
struct index_accessor {
    using element_type     = ElementType;
    using reference        = std::conditional_t<Base == 0,
                                ElementType&, // just like default_accessor when base == 0
                                std::conditional_t<std::is_const_v<ElementType>,
                                    ElementType, // base applied to const accessor directly
                                    index_reference<ElementType,Base>>>; // base applied within index_reference
    using data_handle_type = ElementType*;
    using offset_policy    = index_accessor;

    constexpr index_accessor() = default;

    template<class OtherAccessor, typename = typename std::enable_if_t<
          (!std::is_const_v<ElementType> && std::is_same_v<ElementType, std::remove_const_t<typename OtherAccessor::element_type>>)
        ||( std::is_const_v<ElementType> && std::is_same_v<ElementType, std::add_const_t<typename OtherAccessor::element_type>>)>>
    constexpr index_accessor(const OtherAccessor&) {}

    template<size_t B=Base, typename = std::enable_if_t<B == 0>>
    constexpr ElementType& access(data_handle_type p, size_t i) const noexcept {
        return p[i];
    }
    template<size_t B=Base, typename = std::enable_if_t<B != 0 && !std::is_const_v<ElementType>>>
    constexpr index_reference<ElementType,Base> access(data_handle_type p, size_t i) const noexcept {
        return p[i];
    }
    template<size_t B=Base, typename = std::enable_if_t<B != 0 && std::is_const_v<ElementType>>>
    constexpr ElementType access(data_handle_type p, size_t i) const noexcept {
        return p[i] - base_;
    }

    constexpr data_handle_type offset(data_handle_type p, size_t i) const noexcept {
        return p + i;
    }

    constexpr operator default_accessor<element_type>() const noexcept {
        return default_accessor<element_type>();
    }

    static constexpr ElementType base_{Base};
};

} // namespace atlas

// ------------------------------------------------------------------------------------------------
