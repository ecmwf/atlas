#include "atlas/grid/detail/grid/CubedSphere2.h"

#include <cmath>

#include "eckit/geometry/Sphere.h"
#include "eckit/utils/Hash.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

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

// Get the lonlat for a given index
void CubedSphere2::xy(idx_t n, Point2& point) const {
    auto [t, i, j] = get_cs_indices(n);

    PointAlphaBeta ab = ij_to_curvilinear_coord(i, j);
    PointXY tangent_xy = curvilinear_to_tangent_coord(ab);
    PointXYZ xyz = tangent_to_xyz_coord(tangent_xy, t);
    Matrix xyz_m { {xyz[0], xyz[1], xyz[2]} };

    eckit::geometry::Sphere::convertCartesianToSpherical(xyz_m.norm(), xyz, point);
}

Point2 CubedSphere2::xy(idx_t n) const {
    Point2 point;
    xy(n, point);
    return point;
}

void CubedSphere2::lonlat(idx_t n, Point2& point) const {
    xy(n, point);
    projection_.xy2lonlat(point);
}

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

CubedSphere2::PointAlphaBeta CubedSphere2::ij_to_curvilinear_coord(idx_t i, idx_t j) const {
    const auto get_coord = [&](idx_t idx) {
        return M_PI / 2 * (-0.5 + (0.5 + idx) / N());
    };
    return {get_coord(i), get_coord(j)};
}

PointXY CubedSphere2::curvilinear_to_tangent_coord(CubedSphere2::PointAlphaBeta& curvi_coord) const {
    return {std::tan(curvi_coord[0]), std::tan(curvi_coord[1])};
}

PointXYZ CubedSphere2::tangent_to_xyz_coord(PointXY& tan_coord, idx_t tile) const {
    Matrix tan_point { {tan_coord[0]}, {tan_coord[1]}, {1} };
    
    Matrix xyz(3, 1);
    xyz = lfric_rotations_[tile].transpose() * tan_point;

    return {xyz(0), xyz(1), xyz(2)};
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
