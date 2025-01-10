#pragma once

#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/grid/CubedSphere2.h"

namespace atlas {

class CubedSphereGrid2 : public atlas::Grid {
public:
    using grid_t = grid::detail::grid::CubedSphere2;

public:
    CubedSphereGrid2(idx_t resolution);
    CubedSphereGrid2(const Grid& grid);
    CubedSphereGrid2(idx_t resolution, Projection projection);

    bool valid() const { return grid_; };
    operator bool() const { return valid(); }

private:
    const grid_t* grid_;
};

}  // namespace atlas
