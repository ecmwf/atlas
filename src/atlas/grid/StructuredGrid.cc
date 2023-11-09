/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/Grid.h"

#include <limits>
#include <string>
#include <vector>

#include "eckit/config/Parametrisation.h"

#include "atlas/domain/Domain.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Spacing.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/grid/detail/grid/Gaussian.h"
#include "atlas/projection/Projection.h"
#include "atlas/util/Config.h"

namespace atlas {

inline const StructuredGrid::grid_t* structured_grid(const Grid::Implementation* grid) {
    return dynamic_cast<const StructuredGrid::grid_t*>(grid);
}

StructuredGrid::StructuredGrid(): Grid(), grid_(nullptr) {}

StructuredGrid::StructuredGrid(const Grid& grid): Grid(grid), grid_(structured_grid(get())) {}

StructuredGrid::StructuredGrid(const Grid::Implementation* grid): Grid(grid), grid_(structured_grid(get())) {}

StructuredGrid::StructuredGrid(const std::string& grid, const Domain& domain):
    Grid(grid, domain), grid_(structured_grid(get())) {}

StructuredGrid::StructuredGrid(const std::string& grid, const Projection& projection, const Domain& domain):
    Grid(grid, projection, domain), grid_(structured_grid(get())) {}

StructuredGrid::StructuredGrid(const Config& grid): Grid(grid), grid_(structured_grid(get())) {}

StructuredGrid::StructuredGrid(const XSpace& xspace, const YSpace& yspace, const Projection& projection,
                               const Domain& domain):
    Grid(new StructuredGrid::grid_t(xspace, yspace, projection, domain)), grid_(structured_grid(get())) {}

StructuredGrid::StructuredGrid(const Grid& grid, const Grid::Domain& domain):
    Grid(grid, domain), grid_(structured_grid(get())) {}

ReducedGaussianGrid::ReducedGaussianGrid(const std::vector<long>& nx, const Domain& domain):
    ReducedGaussianGrid::grid_t(grid::detail::grid::reduced_gaussian(nx, domain)) {}

ReducedGaussianGrid::ReducedGaussianGrid(const std::vector<int>& nx, const Domain& domain):
    ReducedGaussianGrid::grid_t(grid::detail::grid::reduced_gaussian(nx, domain)) {}

ReducedGaussianGrid::ReducedGaussianGrid(const std::vector<long>& nx, const Projection& projection):
    ReducedGaussianGrid::grid_t(grid::detail::grid::reduced_gaussian(nx, projection)) {}

ReducedGaussianGrid::ReducedGaussianGrid(const std::vector<int>& nx, const Projection& projection):
    ReducedGaussianGrid::grid_t(grid::detail::grid::reduced_gaussian(nx, projection)) {}

ReducedGaussianGrid::ReducedGaussianGrid(const std::initializer_list<idx_t>& nx):
    ReducedGaussianGrid(std::vector<idx_t>(nx)) {}

RegularGaussianGrid::RegularGaussianGrid(int N, const Grid::Domain& domain):
    RegularGaussianGrid::grid_t("F" + std::to_string(N), domain) {}


inline const HealpixGrid::grid_t* healpix_grid(const Grid::Implementation* grid) {
    return dynamic_cast<const HealpixGrid::grid_t*>(grid);
}

HealpixGrid::HealpixGrid(const Grid& grid): StructuredGrid(grid), grid_(healpix_grid(get())) {}

HealpixGrid::HealpixGrid(int N, const std::string& ordering): HealpixGrid(Grid(new HealpixGrid::grid_t(N, ordering))) {}

}  // namespace atlas
