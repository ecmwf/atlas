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
#include <iostream>
#include <string>

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

namespace projection {
class Jacobian;
}

namespace grid {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace detail {
class CubedSphereTiles;
}  // namespace detail
#endif

//---------------------------------------------------------------------------------------------------------------------

class CubedSphereTiles : DOXYGEN_HIDE(public util::ObjectHandle<atlas::grid::detail::CubedSphereTiles>) {
public:
    using Spec = util::Config;

public:
    using Handle::Handle;
    CubedSphereTiles() = default;
    CubedSphereTiles(const eckit::Parametrisation&);
    CubedSphereTiles(const std::string&);

    /// Type of the cubed-sphere tiles:
    std::string type() const;

    // These are offsets needed for transforming
    // from xy space to the "archetypal base" tile.
    std::array<std::array<double, 6>, 2> xy2abOffsets() const;

    std::array<std::array<double, 6>, 2> ab2xyOffsets() const;

    void rotate(idx_t t, double xyz[]) const;

    void unrotate(idx_t t, double xyz[]) const;

    // tile index from xy space
    idx_t indexFromXY(const double xy[]) const;
    idx_t indexFromXY(const PointXY& xy) const;

    // tile index from longitude and latitude space
    idx_t indexFromLonLat(const double lonlat[]) const;
    idx_t indexFromLonLat(const PointLonLat& lonlat) const;

    // enforceXYdomain reinforces the tile shape in xy space;
    // if values move a miniscule amount outside the domain, it will be brought back in.
    void enforceXYdomain(double xy[]) const;

    idx_t size() const;

    // this provides periodicity to each of the tiles by extending each tile over edges
    // in a cross-like fashion. Periodicity of this form does not allow
    // a "diagonal" extension over corners of the cube.
    atlas::PointXY tileCubePeriodicity(const atlas::PointXY& xyExtended, const atlas::idx_t tile) const;

    /// @brief Return the position of the tile centre in xy space.
    const PointXY& tileCentre(size_t t) const;

    /// @brief Return the Jacobian of xy with respect to the curvilinear
    ///        coordinates of the tile.
    const projection::Jacobian& tileJacobian(size_t t) const;

private:
    /// Output to stream
    void print(std::ostream&) const;

    friend std::ostream& operator<<(std::ostream& s, const CubedSphereTiles& cst);
};

}  // namespace grid
}  // namespace atlas
