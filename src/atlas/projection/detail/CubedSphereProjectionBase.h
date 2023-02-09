/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/grid/Tiles.h"
#include "atlas/library/config.h"
#include "atlas/projection/Jacobian.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "eckit/config/Parametrisation.h"
#include "eckit/utils/Hash.h"

namespace atlas {
class CubedSphereTiles;
}

namespace atlas {
namespace projection {
namespace detail {

class CubedSphereProjectionBase : public ProjectionImpl {
public:
    using ProjectionImpl::jacobian;

    // constructor
    CubedSphereProjectionBase(const eckit::Parametrisation&);

    void hash(eckit::Hash&) const;

    const atlas::grid::CubedSphereTiles& getCubedSphereTiles() const { return tiles_; };

    /// @brief   Convert (x, y) coordinate to (alpha, beta) on tile t.
    ///
    /// @details Converts the Atlas xy coordinates to the angular coordinates
    ///          described of tile t, described by by Ronchi et al. (1996,
    ///          Journal of Computational Physics, 124, 93).
    ///          Note that the xy coordinate must lie within the domain
    ///          (x <= xmin && x >= xmax) || (y <= ymin && y >= ymax) where
    ///          xmin, xmax, ymin and ymax are the boundaries of tile t.
    ///@{
    Point2 xy2alphabeta(const Point2& xy, idx_t t) const;
    virtual void xy2alphabeta(double crd[], idx_t t) const = 0;
    ///@}

    /// @brief   Convert (alpha, beta) coordinate to (x, y) on tile t.
    ///
    /// @details Performs the inverse of xy2alpha beta. Note that the result is
    ///          degenerate when abs(alpha) > 45 && abs(beta) > 45 &&
    ///          abs(alpha) == abs(beta). In these circumstances, the method
    ///          will return the anticlockwise-most of the two possible
    ///          values.
    ///@{
    Point2 alphabeta2xy(const Point2& alphabeta, idx_t t) const;
    virtual void alphabeta2xy(double crd[], idx_t) const = 0;
    ///@}

    /// @brief   Convert (lon, lat) coordinate to (alpha, beta) on tile t
    ///
    /// @details Converts the lon lat coordinates to the angular coordinates
    ///          described of tile t, described by by Ronchi et al. (1996,
    ///          Journal of Computational Physics, 124, 93).
    /// @{
    Point2 lonlat2alphabeta(const Point2& lonlat, idx_t t) const;
    virtual void lonlat2alphabeta(double crd[], idx_t) const = 0;
    /// @}

    /// @brief   Convert (lon, lat) coordinate to (alpha, beta) on tile t
    ///
    /// @details Performs the inverse of lonlat2alphabeta.
    ///
    Point2 alphabeta2lonlat(const Point2& alphabeta, idx_t t) const;
    virtual void alphabeta2lonlat(double crd[], idx_t) const = 0;

    /// @brief   Jacobian of (x, y) with respect to (lon, lat) on tile t
    ///
    /// @details Returns the Jacobian
    ///                     ∂x/∂λ, ∂x/∂φ
    ///                     ∂y/∂λ, ∂y/∂φ
    ///          for tile t.
    virtual Jacobian jacobian(const PointLonLat& lonlat, idx_t t) const = 0;

    /// @brief   Jacobian of (alpha, beta) with respect to (lon, lat) on tile t
    ///
    /// @details Returns the Jacobian
    ///                     ∂α/∂λ, ∂α/∂φ
    ///                     ∂β/∂λ, ∂β/∂φ
    ///          for tile t.
    virtual Jacobian alphabetaJacobian(const PointLonLat& lonlat, idx_t t) const = 0;

protected:
    // projection and inverse projection
    void xy2lonlat_post(double xyz[], const idx_t t, double crd[]) const;
    void lonlat2xy_pre(double crd[], idx_t& t, double xyz[]) const;

    void xy2alphabetat(const double xy[], idx_t& t, double ab[]) const;
    void alphabetat2xy(const idx_t t, const double ab[], double xy[]) const;

private:
    atlas::grid::CubedSphereTiles tiles_;
    // Shift entire grid
    double shiftLon_;
    // Schmidt transform
    bool doSchmidt_;
    double stretchFac_;
    double targetLon_;
    double targetLat_;

    std::array<std::array<double, 6>, 2> tiles_offsets_ab2xy_;
    std::array<std::array<double, 6>, 2> tiles_offsets_xy2ab_;
};

}  // namespace detail
}  // namespace projection
}  // namespace atlas
