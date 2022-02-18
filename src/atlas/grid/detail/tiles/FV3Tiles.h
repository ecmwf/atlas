/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>

#include "atlas/grid/detail/tiles/Tiles.h"

namespace atlas {
namespace grid {
namespace detail {

class FV3CubedSphereTiles : public CubedSphereTiles {
public:
    // constructor
    FV3CubedSphereTiles(const eckit::Parametrisation&);

    static std::string static_type() { return "cubedsphere_fv3"; }

    virtual std::string type() const override { return static_type(); }

    virtual std::array<std::array<double, 6>, 2> xy2abOffsets() const override;

    virtual std::array<std::array<double, 6>, 2> ab2xyOffsets() const override;

    virtual void rotate(idx_t t, double xyz[]) const override;

    virtual void unrotate(idx_t, double xyz[]) const override;

    virtual idx_t indexFromXY(const double xy[]) const override;

    virtual idx_t indexFromLonLat(const double lonlat[]) const override;

    virtual void enforceXYdomain(double xy[]) const override;

    virtual atlas::PointXY tileCubePeriodicity(const atlas::PointXY& xyExtended,
                                               const atlas::idx_t tile) const override;

    virtual void print(std::ostream&) const override;

    virtual const PointXY& tileCentre(size_t t) const override;

    virtual const Jacobian& tileJacobian(size_t t) const override;


private:
};

}  // namespace detail
}  // namespace grid
}  // namespace atlas
