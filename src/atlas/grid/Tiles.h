/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include <string>
#include <iostream>

#include "atlas/library/config.h"
#include "atlas/util/ObjectHandle.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
class Hash;
}  // namespace eckit

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
class PointXY;
class PointLonLat;

namespace util {
class Config;
}  // namespace util

namespace grid {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace detail {
class CubedSphereTiles;
}  // namespace detail
#endif

//---------------------------------------------------------------------------------------------------------------------

class CubedSphereTiles : DOXYGEN_HIDE(
    public util::ObjectHandle<atlas::grid::detail::CubedSphereTiles> ) {
public:
    using Spec = util::Config;

public:
    using Handle::Handle;
    CubedSphereTiles() = default;
    CubedSphereTiles( const eckit::Parametrisation& );
    CubedSphereTiles( const std::string& );

    /// Type of the cubed-sphere tiles:
    std::string type() const;

    std::array<std::array<double,6>,2> xy2abOffsets() const;

    std::array<std::array<double,6>,2> ab2xyOffsets() const;

    void rotate( idx_t t, double xyz[] ) const;

    void unrotate( idx_t t, double xyz[] ) const;

    idx_t indexFromXY( const double xy[] ) const;

    idx_t indexFromLonLat( const double lonlat[] ) const;

    void enforceXYdomain( double xy[] ) const;

    idx_t size() const;

private:
    /// Output to stream
    void print( std::ostream& ) const;

    friend std::ostream& operator<<( std::ostream& s, const CubedSphereTiles& cst );

};

}  // namespace grid
}  // namespace atlas
