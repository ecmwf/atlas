/*
 * (C) Crown Copyright 2021, Met Office.
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

#include "atlas/grid/detail/tiles/Tiles.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace grid {
namespace detail {

class LFRicCubedSphereTiles : public CubedSphereTiles {
public:
    // constructor
    LFRicCubedSphereTiles(const eckit::Parametrisation&);

    static std::string static_type() { return "cubedsphere_lfric"; }

    virtual std::string type() const override { return static_type(); }

    virtual std::array<std::array<double, 6>, 2> xy2abOffsets() const override;

    virtual std::array<std::array<double, 6>, 2> ab2xyOffsets() const override;

    virtual void rotate(idx_t t, double xyz[]) const override;

    virtual void unrotate(idx_t t, double xyz[]) const override;

    virtual idx_t indexFromXY(const double xy[]) const override;

    virtual idx_t indexFromLonLat(const double lonlat[]) const override;

    virtual void enforceXYdomain(double xy[]) const override;

    virtual PointXY tileCubePeriodicity(const PointXY& xyExtended, const idx_t tile) const override;

    virtual void print(std::ostream&) const override;

    virtual const PointXY& tileCentre(size_t t) const override;

    virtual const Jacobian& tileJacobian(size_t t) const override;

private:
    static PointXY botLeftTile(size_t t);
    static PointXY botRightTile(size_t t);
    static PointXY topLeftTile(size_t t);
    static PointXY topRightTile(size_t t);

    bool withinCross(const idx_t t, const PointXY& withinRange) const;

    void enforceWrapAround(const idx_t t, PointXY& withinRange) const;


    // Centre of each tile in xy-space.
    static const std::array<PointXY, 6> tileCentres_;
    // Jacobian of xy with respect to tile curvilinear coordinates.
    static const std::array<Jacobian, 6> tileJacobians_;
};

}  // namespace detail
}  // namespace grid
}  // namespace atlas
