#include "atlas/grid/detail/grid/CubedSphere2.h"

#include <cmath>

#include "eckit/geometry/Sphere.h"
#include "eckit/utils/Hash.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

// Public methods

CubedSphere2::CubedSphere2(idx_t resolution) : N_(resolution) {}

std::string CubedSphere2::name() const {
    return "CS-LFR-" + std::to_string(N_) + "-2";
}

std::string CubedSphere2::type() const {
    return type_;
}

// Provide a unique identification hash for the grid and the projection.
void CubedSphere2::hash(eckit::Hash& h) const {
    h.add(name());
    h.add(int(N_));

    // also add projection information
    projection().hash(h);

    // also add domain information, even though already encoded in grid.
    domain().hash(h);
}

// Return the bounding box for the grid, global
RectangularLonLatDomain CubedSphere2::lonlatBoundingBox() const {
    return GlobalDomain();
}

// Return the total number of points
idx_t CubedSphere2::size() const {
    return N_ * N_ * nTiles_;
}

// Return the specification for the grid.
Grid::Spec CubedSphere2::spec() const {
    Grid::Spec grid_spec;

    grid_spec.set("name", name());
    grid_spec.set("type", type());
    grid_spec.set("projection", projection().spec());
    grid_spec.set("domain", domain());

    return grid_spec;
}

// Get the xy for a given index
void CubedSphere2::xy(idx_t n, Point2& point) const {
    auto [t, i, j] = get_cs_indices(n);

    PointXY tangent_xy = ij_to_tangent_coord(i, j);
    PointXYZ xyz = tangent_to_xyz_coord(tangent_xy, t);
    eckit::geometry::Sphere::convertCartesianToSpherical(1., xyz, point);
}

// Get the xy for a given index
Point2 CubedSphere2::xy(idx_t n) const {
    Point2 point;
    xy(n, point);
    return point;
}

// Get the lonlat for a given index
void CubedSphere2::lonlat(idx_t n, Point2& point) const {
    xy(n, point);
    projection_.xy2lonlat(point);
}

// Get the lonlat for a given index
Point2 CubedSphere2::lonlat(idx_t n) const {
    Point2 point;
    lonlat(n, point);
    return point;
}

// Protected methods

// Print the name of the Grid
void CubedSphere2::print(std::ostream& os) const {
    os << "CubedSphere2(Name:" << name() << ")";
}

// Private methods

// Get t, i, and j for a given index
CubedSphere2::CSIndices CubedSphere2::get_cs_indices(gidx_t n) const {
    ATLAS_ASSERT(n <= size());
    const idx_t tile_size = N() * N();
    const idx_t t = n / tile_size;
    const idx_t ij = n % tile_size;
    const idx_t j = ij / N();
    const idx_t i = ij % N();
    return {t, i, j};
}

// Get the point on the tangent plane for a given ij index
PointXY CubedSphere2::ij_to_tangent_coord(idx_t i, idx_t j) const {
    const auto get_curvilinear_coord = [&](idx_t idx) {
        return M_PI_2 * (-0.5 + (0.5 + static_cast<double>(idx)) / static_cast<double>(N()));
    };
    return {std::tan(get_curvilinear_coord(i)), std::tan(get_curvilinear_coord(j))};
}

// Transform a point on the tangent plane to a point on a cube
PointXYZ CubedSphere2::tangent_to_xyz_coord(const PointXY& tan_coord, idx_t tile) const {
    PointXYZ xyz;
    const Matrix& transform = lfric_rotations_transposed_[tile];

    xyz[0] = transform[0][0] * tan_coord[0] + transform[0][1] * tan_coord[1] + transform[0][2];
    xyz[1] = transform[1][0] * tan_coord[0] + transform[1][1] * tan_coord[1] + transform[1][2];
    xyz[2] = transform[2][0] * tan_coord[0] + transform[2][1] * tan_coord[1] + transform[2][2];

    return PointXYZ::normalize(xyz);
}

std::string CubedSphere2::name() const {
    return "CS-LFR-" + std::to_string(N_) + "-2";
}

std::string CubedSphere2::type() const {
    return type_;
}

// Provide a unique identification hash for the grid and the projection.
void CubedSphere2::hash(eckit::Hash& h) const {
    h.add(name()); // use name() or type()?
    h.add(int(N_));

    // also add projection information
    projection().hash(h);

    // also add domain information, even though already encoded in grid.
    domain().hash(h);
}

// Return the bounding box for the grid, global
RectangularLonLatDomain CubedSphere2::lonlatBoundingBox() const {
    return GlobalDomain();
}

idx_t CubedSphere2::size() const {
    // Example from CubedSphere.h
    // return accumulate(npts_.begin(), npts_.end(), 0);

    // For now, return expected size
    return N_ * N_ * nTiles_;
}

// Return the specification for the grid.
Grid::Spec CubedSphere2::spec() const {
    // Copied from CubedSphere.cc
    Grid::Spec grid_spec;

    if (type() == "cubedsphere2") {
        grid_spec.set("name", name());
    }
    else {
        grid_spec.set("type", type());
    }
    grid_spec.set("projection", projection().spec());
    return grid_spec;
}

// Print the name of the Grid
void CubedSphere2::print(std::ostream& os) const {
    os << "CubedSphere2(Name:" << name() << ")";
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
