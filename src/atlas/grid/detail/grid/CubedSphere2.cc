#include "atlas/grid/detail/grid/CubedSphere2.h"

#include <cmath>

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

    // 1. Get point on base face (xy plane)
    Matrix base_point(3, 1);
    base_point(0) = std::tan(index_to_curvilinear(i));
    base_point(1) = std::tan(index_to_curvilinear(j));
    base_point(2) = 1;

    // 2. Apply rotation (move point from xy plane to 3D cube)
    Matrix xyz(3, 1);
    xyz = lfric_rotations_[t].transpose() * base_point;

    // // 3. Project the point onto a (cubed)sphere
    point[0] = std::atan2(xyz(1), xyz(0)) * rad_to_deg_;
    point[1] = std::asin(xyz(2) / (
        std::sqrt(xyz(0)*xyz(0) + xyz(1)*xyz(1) + xyz(2)*xyz(2)) // Magnitude
        )) * rad_to_deg_;
}

Point2 CubedSphere2::xy(idx_t n) const {
    Point2 point;
    xy(n, point);
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

// Get point between [-pi/4..pi/4] given index and number of points along edge
// Applies offset to get from prime to dual mesh
double CubedSphere2::index_to_curvilinear(idx_t index) const {
    return M_PI / 2 * (-0.5 + (0.5 + index) / N_);
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
