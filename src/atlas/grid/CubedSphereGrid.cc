/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <string>

#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/grid/CubedSphere.h"
#include "atlas/projection/Projection.h"

namespace atlas {

inline const CubedSphereGrid::grid_t* cubedsphere_grid(const Grid::Implementation* grid) {
    return dynamic_cast<const CubedSphereGrid::grid_t*>(grid);
}

CubedSphereGrid::CubedSphereGrid(): Grid(), grid_(nullptr) {}

CubedSphereGrid::CubedSphereGrid(const Grid& grid): Grid(grid), grid_(cubedsphere_grid(get())) {}

CubedSphereGrid::CubedSphereGrid(const Grid::Implementation* grid): Grid(grid), grid_(cubedsphere_grid(get())) {}

CubedSphereGrid::CubedSphereGrid(const std::string& grid): Grid(grid), grid_(cubedsphere_grid(get())) {}

CubedSphereGrid::CubedSphereGrid(const int& N, const Projection& projection):
    Grid(new CubedSphereGrid::grid_t(N, projection, "C")), grid_(cubedsphere_grid(get())) {}

}  // namespace atlas
