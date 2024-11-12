#include "atlas/grid/CubedSphereGrid2.h"
#include "atlas/grid/detail/grid/CubedSphere2.h"

namespace atlas {

CubedSphereGrid2::CubedSphereGrid2(idx_t resolution)
    : Grid(new grid::detail::grid::CubedSphere2(resolution)) {}

}  // namespace atlas
