/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/library/config.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Point.h"

namespace atlas {
class CubedSphereGrid;
}

namespace atlas {}

namespace atlas {
namespace projection {
namespace detail {
class CubedSphereProjectionBase;
}  // namespace detail
}  // namespace projection
}  // namespace atlas

namespace atlas {
namespace meshgenerator {
namespace detail {
namespace cubedsphere {

/// Enum for (i, j, t) coordinate fields.
struct Coordinates {
    enum k : idx_t
    {
        T,
        I,
        J
    };
};

/// Enum for tile edges.
struct TileEdge {
    enum k : size_t
    {
        LEFT,
        BOTTOM,
        RIGHT,
        TOP,
        UNDEFINED
    };
};

/// Class to store (i, j) indices as a Point2 coordinate.
class PointIJ : public Point2 {
public:
    using Point2::Point2;
    PointIJ(): Point2() {}

    /// Index constructor.
    template <typename IndexI, typename IndexJ>
    inline PointIJ(IndexI i, IndexJ j): Point2(static_cast<double>(i), static_cast<double>(j)) {}

    /// @{
    ///  Return i or j by value.
    inline double i() const { return x_[0]; }
    inline double j() const { return x_[1]; }
    /// @}

    /// @{
    /// Return i or j by reference
    inline double& i() { return x_[0]; }
    inline double& j() { return x_[1]; }
    /// @}

    /// @{
    /// Round i or j to node index.
    inline idx_t iNode() const { return static_cast<idx_t>(std::round(i())); }
    inline idx_t jNode() const { return static_cast<idx_t>(std::round(j())); }
    /// @}

    /// @{
    /// Round i or j to cell index.
    inline idx_t iCell() const { return static_cast<idx_t>(std::floor(i())); }
    inline idx_t jCell() const { return static_cast<idx_t>(std::floor(j())); }
    /// @}
};

/// (t, PointXY) tuple.
class PointTXY : public std::pair<idx_t, PointXY> {
    using std::pair<idx_t, PointXY>::pair;

public:
    idx_t& t() { return first; }
    PointXY& xy() { return second; }
    const idx_t& t() const { return first; }
    const PointXY& xy() const { return second; }
};

/// (t, PointIJ) tuple.
class PointTIJ : public std::pair<idx_t, PointIJ> {
    using std::pair<idx_t, PointIJ>::pair;

public:
    idx_t& t() { return first; }
    PointIJ& ij() { return second; }
    const idx_t& t() const { return first; }
    const PointIJ& ij() const { return second; }
};

/// \brief   Class to convert between ij and xy on a tile and its four
///          surrounding neighbours.
///
///          Helper class to deal with the coordinate system roations and
///          displacements between a tile and its neighbours. This class
///          is specifcially written to comupute the (x, y) and (t, i, j)
///          coordinates of halos that extend across tile boundaries.
class NeighbourJacobian {
private:
    using Jacobian = projection::Jacobian;

public:
    /// Default constructor.
    NeighbourJacobian() = default;

    /// Grid-data constructor.
    NeighbourJacobian(const CubedSphereGrid& csGrid);

    /// Convert ij on local tile t to xy.
    PointXY xy(const PointIJ& ij, idx_t t) const;

    /// Convert xy on local tile t to ij.
    PointIJ ij(const PointXY& xy, idx_t t) const;

    /// Convert extrapolated xy on tile t to global xy and t (needed for halos).
    PointTXY xyLocalToGlobal(const PointXY& xyLocal, idx_t tLocal) const;

    /// Convert extrapolated ij on tile t to global ij and t (needed for halos).
    PointTIJ ijLocalToGlobal(const PointIJ& ijLocal, idx_t tLocal) const;

    /// Return true if ij is interior or on the edge of a tile.
    bool ijInterior(const PointIJ& ij) const;

    /// Return true if ij is on the edge of a tile.
    bool ijEdge(const PointIJ& ij) const;

    /// Return true if ij is in the valid "+" halo extension of at tile.
    bool ijCross(const PointIJ& ij) const;

    /// Makes sure points near tile edges are *exactly* on the edge.
    PointXY snapToEdge(const PointXY& xy, idx_t t) const;

private:
    // Pointer to grid projection.
    const projection::detail::CubedSphereProjectionBase* csProjection_{};

    // Grid size.
    idx_t N_{};

    // Jacobian of xy with respect to ij for each tile.
    std::array<Jacobian, 6> dxy_by_dij_{};

    // Jacobian of ij with respect to xy for each tile.
    std::array<Jacobian, 6> dij_by_dxy_{};

    // Lower-left xy position on each tile.
    std::array<PointXY, 6> xy00_{};

    // Min xy on each tile.
    std::array<PointXY, 6> xyMin_{};

    // Max xy on each tile.
    std::array<PointXY, 6> xyMax_{};

    // Properties of four neighbours of a tile.
    struct Neighbours {
        // Tile ID.
        std::array<idx_t, 4> t_{};

        // Jacobian of remote xy with respect to local xy.
        std::array<Jacobian, 4> dxyGlobal_by_dxyLocal_{};

        // Lower left most local xy position on neighbour tiles.
        std::array<PointXY, 4> xy00Local_;
        std::array<PointXY, 4> xy00Global_;
    };

    // Set of neighbours for each tile.
    std::array<Neighbours, 6> neighbours_{};
};

}  // namespace cubedsphere
}  // namespace detail
}  // namespace meshgenerator
}  // namespace atlas
