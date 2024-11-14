#include "atlas/grid/detail/grid/CubedSphere2.h"

#include <cmath>

#include "eckit/utils/Hash.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

CubedSphere2::CubedSphere2(idx_t resolution) : N_(resolution) {}

std::string CubedSphere2::name() const {
    // - - - - - - TEST - - - - - -
    // Temporarily here to test lonlat()
    PointLonLat point;
    for (int i = 0; i < size(); ++i) {
        lonlat(i, point);
        std::cout << "[" << point[0] << ", " << point[1] << "], ";
    }
    std::cout << std::endl;
    // - - - - - END TEST - - - - -

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

// Return the total number of points
idx_t CubedSphere2::size() const {
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

// Get the lonlat for a given index
void CubedSphere2::lonlat(idx_t n, PointLonLat& point) const {
    // 1. Get point on base face (xy plane)
    int ij = get_tij(n);
    double base_point[3];
    base_point[0] = std::tan(index_to_curvilinear(get_ti(ij)));
    base_point[1] = std::tan(index_to_curvilinear(get_tj(ij)));
    base_point[2] = 1;

    // 2. Apply rotation (move point from xy plane to 3D cube)
    int tile = get_tile(n);
    double xyz[3];
    xyz[0] =  lfric_rotations_[tile][0] * base_point[0]
            + lfric_rotations_[tile][1] * base_point[1]
            + lfric_rotations_[tile][2] * base_point[2];
    xyz[1] =  lfric_rotations_[tile][3] * base_point[0]
            + lfric_rotations_[tile][4] * base_point[1]
            + lfric_rotations_[tile][5] * base_point[2];
    xyz[2] =  lfric_rotations_[tile][6] * base_point[0]
            + lfric_rotations_[tile][7] * base_point[1]
            + lfric_rotations_[tile][8] * base_point[2];

    // 3. Project the point onto a (cubed)sphere
    point[0] = std::atan2(xyz[1], xyz[0]) * rad_to_deg_;
    point[1] = std::asin(xyz[2] / (
        std::sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]) // Magnitude
        )) * rad_to_deg_;
}

PointLonLat CubedSphere2::lonlat(idx_t n) const {
    PointLonLat point;
    lonlat(n, point);
    return point;
}

// Protected methods

// Print the name of the Grid
void CubedSphere2::print(std::ostream& os) const {
    os << "CubedSphere2(Name:" << name() << ")";
}

// Private methods

// Get tile of point given global index
int CubedSphere2::get_tile(idx_t n) const {
    if (n >= N_ * N_ * 6) return -1;
    return n / (N_ * N_);
}

// Get index of point on tile given global index
int CubedSphere2::get_tij(idx_t n) const {
    return n % (N_ * N_);
}

// Get column (lon) of point on tile given global index
int CubedSphere2::get_ti(idx_t n) const {
    return get_tij(n) % N_;
}

// Get row (lat) of point on tile given global index
int CubedSphere2::get_tj(idx_t n) const {
    return get_tij(n) / N_;
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