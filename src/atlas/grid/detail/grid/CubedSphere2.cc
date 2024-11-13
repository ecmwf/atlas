#include "atlas/grid/detail/grid/CubedSphere2.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

CubedSphere2::CubedSphere2(idx_t resolution)
{

}

std::string CubedSphere2::name() const {
    return "CS-LFR-" + std::to_string(N_) + "-2";
}

std::string CubedSphere2::type() const {
    return type_;
}

// Provide a unique identification hash for the grid and the projection.
void CubedSphere2::hash(eckit::Hash& h) const {
    h.add("CubedSphere2");
    h.add(int(N_));

    // also add projection information
    projection().hash(h);

    // also add domain information, even though already encoded in grid.
    domain().hash(h);
}

// Return the bounding box for the grid, global
RectangularLonLatDomain CubedSphere2::lonlatBoundingBox() const {
    return projection_ ? projection_.lonlatBoundingBox(computeDomain()) : domain();
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

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas