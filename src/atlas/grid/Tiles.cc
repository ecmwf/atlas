/*
 * (C) Crown Copyright 2021 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/Tiles.h"
#include "atlas/grid/detail/tiles/FV3Tiles.h"
#include "atlas/grid/detail/tiles/LFRicTiles.h"
#include "atlas/grid/detail/tiles/Tiles.h"


namespace atlas {
namespace grid {

CubedSphereTiles::CubedSphereTiles(const eckit::Parametrisation& p):
    Handle(atlas::grid::detail::CubedSphereTiles::create(p)) {}

CubedSphereTiles::CubedSphereTiles(const std::string& s): Handle(atlas::grid::detail::CubedSphereTiles::create(s)) {}


std::string CubedSphereTiles::type() const {
    return get()->type();
}

std::array<std::array<double, 6>, 2> CubedSphereTiles::xy2abOffsets() const {
    return get()->xy2abOffsets();
}

std::array<std::array<double, 6>, 2> CubedSphereTiles::ab2xyOffsets() const {
    return get()->ab2xyOffsets();
}

void CubedSphereTiles::rotate(idx_t t, double xyz[]) const {
    return get()->rotate(t, xyz);
}

void CubedSphereTiles::unrotate(idx_t t, double xyz[]) const {
    return get()->unrotate(t, xyz);
}

idx_t CubedSphereTiles::indexFromXY(const double xy[]) const {
    return get()->indexFromXY(xy);
}

idx_t CubedSphereTiles::indexFromLonLat(const double lonlat[]) const {
    return get()->indexFromLonLat(lonlat);
}

idx_t CubedSphereTiles::size() const {
    return get()->size();
}

void CubedSphereTiles::enforceXYdomain(double xy[]) const {
    return get()->enforceXYdomain(xy);
}

atlas::PointXY CubedSphereTiles::tileCubePeriodicity(const atlas::PointXY& xyExtended, const atlas::idx_t tile) const {
    return get()->tileCubePeriodicity(xyExtended, tile);
}

void CubedSphereTiles::print(std::ostream& os) const {
    get()->print(os);
}

std::ostream& operator<<(std::ostream& os, const CubedSphereTiles& t) {
    t.print(os);
    return os;
}

const PointXY& CubedSphereTiles::tileCentre(size_t t) const {
    return get()->tileCentre(t);
}

const Jacobian& CubedSphereTiles::tileJacobian(size_t t) const {
    return get()->tileJacobian(t);
}

}  // namespace grid
}  // namespace atlas
