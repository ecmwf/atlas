/*
 * (C) Copyright 2026- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

 #pragma once

#include <cstdint>
#include <type_traits>

namespace atlas {

class GridPointFlags {
public:
    using Representation = int;

    enum Bit : Representation {
        none      = 0,
        duplicate = 1 << 0,
        extension = 1 << 1,
        invalid   = 1 << 2
    };

    constexpr GridPointFlags() = default;
    constexpr GridPointFlags(Bit flag): value_(flag) {}
    explicit constexpr GridPointFlags(Representation value): value_(value) {}

    constexpr Representation value() const { return value_; }
    constexpr bool has(GridPointFlags flags) const { return (value_ & flags.value_) == flags.value_; }

    constexpr GridPointFlags& set(GridPointFlags flags) {
        value_ |= flags.value_;
        return *this;
    }

    constexpr GridPointFlags& unset(GridPointFlags flags) {
        value_ &= ~flags.value_;
        return *this;
    }

    constexpr GridPointFlags& operator|=(GridPointFlags flags) { return set(flags); }

    friend constexpr GridPointFlags operator|(GridPointFlags lhs, GridPointFlags rhs) { return lhs.set(rhs); }
    friend constexpr bool operator==(GridPointFlags lhs, GridPointFlags rhs) { return lhs.value_ == rhs.value_; }
    friend constexpr bool operator!=(GridPointFlags lhs, GridPointFlags rhs) { return !(lhs == rhs); }

private:
    Representation value_{0};
};

static_assert(std::is_standard_layout<GridPointFlags>::value, "GridPointFlags must have a stable scalar layout");
static_assert(std::is_trivially_copyable<GridPointFlags>::value, "GridPointFlags must be transferable as bytes");
static_assert(sizeof(GridPointFlags) == sizeof(GridPointFlags::Representation), "GridPointFlags must have the same representation as int");
static_assert(alignof(GridPointFlags) == alignof(GridPointFlags::Representation), "GridPointFlags must have the same alignment as int");

}  // namespace atlas