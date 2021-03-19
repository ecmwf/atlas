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

#include <cstdint>

#include "eckit/types/SemanticVersion.h"

namespace atlas {
namespace io {

struct Version {             // 16 bytes
    std::uint32_t major{0};  ///<  Major version
    std::uint32_t minor{1};  ///<  Minor version
    std::uint32_t patch{0};  ///<  Patch version
    std::uint32_t tweak{0};  ///<  Tweak version

    std::string str() const {
        return std::to_string( major ) + "." + std::to_string( minor ) + "." + std::to_string( patch );
    }
    operator std::string() const { return str(); }
    operator eckit::SemanticVersion() const { return eckit::SemanticVersion{major, minor, patch}; }

    bool operator<( const Version& v ) const {
        return eckit::SemanticVersion{major, minor, patch} < eckit::SemanticVersion{v.major, v.minor, v.patch};
    }
    bool operator==( const Version& v ) const {
        return eckit::SemanticVersion{major, minor, patch} == eckit::SemanticVersion{v.major, v.minor, v.patch};
    }
    bool operator!=( const Version& v ) const { return !( *this == v ); }
    bool operator<=( const Version& v ) const { return ( *this < v ) or ( *this == v ); }
    bool operator>( const Version& v ) const { return !( *this <= v ); }
    bool operator>=( const Version& v ) const { return ( *this > v ) or ( *this == v ); }


    friend std::ostream& operator<<( std::ostream& out, const Version& v ) {
        out << v.str();
        return out;
    }
};


}  // namespace io
}  // namespace atlas
