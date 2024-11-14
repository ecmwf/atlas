#include "atlas/grid/CubedSphereGrid2.h"
#include "atlas/grid/detail/grid/CubedSphere2.h"

namespace atlas {

inline const CubedSphereGrid2::grid_t* cubedsphere_grid2(const Grid::Implementation* grid) {
    return dynamic_cast<const CubedSphereGrid2::grid_t*>(grid);
}

CubedSphereGrid2::CubedSphereGrid2(idx_t resolution):
    Grid(new grid::detail::grid::CubedSphere2(resolution)), grid_(cubedsphere_grid2(get())) {}

CubedSphereGrid2::CubedSphereGrid2(const Grid& grid):
    Grid(grid), grid_(cubedsphere_grid2(get())) {}

// Temporarily here for testing lonlat()
void CubedSphereGrid2::printCubedSphere2() const {
    get()->name();
}

}  // namespace atlas
