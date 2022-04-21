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

#include "atlas/library/config.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/util/Config.h"
#include "atlas/util/Object.h"

namespace eckit {
class Parametrisation;
}


namespace atlas {
class PointXY;
}  // namespace atlas

namespace atlas {
namespace grid {
namespace detail {

class CubedSphereTiles : public util::Object {
public:
    using Spec     = util::Config;
    using Jacobian = projection::Jacobian;

public:
    static const CubedSphereTiles* create();

    static const CubedSphereTiles* create(const eckit::Parametrisation&);

    static const CubedSphereTiles* create(const std::string&);

    virtual std::string type() const = 0;

    virtual std::array<std::array<double, 6>, 2> xy2abOffsets() const = 0;

    virtual std::array<std::array<double, 6>, 2> ab2xyOffsets() const = 0;

    virtual void rotate(idx_t t, double xyz[]) const = 0;

    virtual void unrotate(idx_t t, double xyz[]) const = 0;

    virtual idx_t indexFromXY(const double xy[]) const = 0;

    virtual idx_t indexFromLonLat(const double lonlat[]) const = 0;

    idx_t indexFromXY(const PointXY& xy) const;

    idx_t indexFromLonLat(const PointLonLat& lonlat) const;

    virtual void enforceXYdomain(double xy[]) const = 0;

    virtual const PointXY& tileCentre(size_t t) const = 0;

    virtual const Jacobian& tileJacobian(size_t t) const = 0;

    idx_t size() const { return 6; }

    virtual atlas::PointXY tileCubePeriodicity(const atlas::PointXY& xyExtended, const atlas::idx_t tile) const = 0;

    /// Output to stream
    virtual void print(std::ostream&) const = 0;

    friend std::ostream& operator<<(std::ostream& s, const CubedSphereTiles& cst) {
        cst.print(s);
        return s;
    }
};

}  // namespace detail
}  // namespace grid
}  // namespace atlas
