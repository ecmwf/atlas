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

#include <string>

namespace atlas {
namespace util {

//----------------------------------------------------------------------------------------------------------------------

/// @brief Convenience class to modify and interpret bitflags
class Bitflags {
public:
    static void reset(int& flags, int bit = 0) { flags = bit; }

    static void set(int& flags, int bit) { flags |= bit; }

    static void unset(int& flags, int bit) { flags &= (~bit); }

    static void toggle(int& flags, int bit) { flags ^= bit; }

    static bool check(int flags, int bits) { return (flags & bits) == bits; }

    static bool check_all(int flags, int bits) { return (flags & bits) == bits; }

    static bool check_any(int flags, int bits) { return flags & bits; }

    static std::string bitstr(int flags) {
        char str[9] = {0};
        int i;
        for (i = 7; i >= 0; i--) {
            str[i] = (flags & 1) ? '1' : '0';
            flags >>= 1;
        }
        return std::string(str, 9);
    }

    /// @brief Create convenience accessor to modify flags
    /// @note `auto` for return type! BitflagsView is a detail
    static auto view(int& flags);

    /// @brief Create convenience accessor to modify flags
    /// @note `auto` for return type! BitflagsView is a detail
    static auto view(const int& flags);
};

//----------------------------------------------------------------------------------------------------------------------

namespace detail {

//----------------------------------------------------------------------------------------------------------------------

template <typename T>
class BitflagsView {};

//----------------------------------------------------------------------------------------------------------------------

// Template specialication for constant flags. There are no functions to edit the flags
template <>
class BitflagsView<const int> {
    const int flags_;

public:
    BitflagsView(const int flags): flags_(flags) {}
    bool check(int bit) const { return Bitflags::check(flags_, bit); }
    bool check_all(int bit) const { return Bitflags::check_all(flags_, bit); }
    bool check_any(int bit) const { return Bitflags::check_any(flags_, bit); }
};

//----------------------------------------------------------------------------------------------------------------------

// Template specialication for nonconst flags. There are functions to edit the flags
template <>
class BitflagsView<int> {
    int& flags_;

public:
    BitflagsView(int& flags): flags_(flags) {}
    void reset(int bit = 0) { Bitflags::reset(flags_, bit); }
    void set(int bit) { Bitflags::set(flags_, bit); }
    void unset(int bit) { Bitflags::unset(flags_, bit); }
    void toggle(int bit) { Bitflags::toggle(flags_, bit); }
    bool check(int bit) const { return Bitflags::check(flags_, bit); }
    bool check_all(int bit) const { return Bitflags::check_all(flags_, bit); }
    bool check_any(int bit) const { return Bitflags::check_any(flags_, bit); }
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace detail

//----------------------------------------------------------------------------------------------------------------------

inline auto Bitflags::view(const int& flags) {
    return detail::BitflagsView<const int>(flags);
}
inline auto Bitflags::view(int& flags) {
    return detail::BitflagsView<int>(flags);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
