/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grid/Grid.h"

#include <limits>
#include <vector>
#include "eckit/config/Parametrisation.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Spacing.h"
#include "atlas/grid/Domain.h"
#include "atlas/projection/Projection.h"
#include "atlas/util/Config.h"
#include "atlas/grid/detail/grid/Gaussian.h"
#include "atlas/grid/detail/grid/Structured.h"

namespace atlas {

Grid::Grid():
    grid_( nullptr ) {
}

Grid::Grid(const Grid& grid):
    grid_( grid.grid_ ) {
}

Grid::Grid( const Grid::Implementation *grid ):
    grid_( grid ) {
}

Grid::Grid( const std::string& shortname, const Domain& domain ) {
    grid_ = Grid::Implementation::create(
                shortname, 
                Config("domain", domain.spec())
            );
}

Grid::Grid( const Config& p ) {
    grid_ = Grid::Implementation::create(p);
}

namespace grid {

inline const UnstructuredGrid::grid_t* unstructured_grid( const Grid::Implementation *grid ) {
    return dynamic_cast<const UnstructuredGrid::grid_t*>(grid);
}

UnstructuredGrid::UnstructuredGrid():
    Grid(),
    grid_( nullptr ) {
}

UnstructuredGrid::UnstructuredGrid( const Grid& grid ):
    Grid( grid ),
    grid_( unstructured_grid(get()) ) {
}

UnstructuredGrid::UnstructuredGrid( const Grid::Implementation* grid ):
    Grid( grid ),
    grid_( unstructured_grid(get()) ) {
}

UnstructuredGrid::UnstructuredGrid( const Config& grid ):
    Grid( grid ),
    grid_( unstructured_grid(get()) ) {
}

UnstructuredGrid::UnstructuredGrid( std::vector<PointXY>* xy ):
    Grid( new UnstructuredGrid::grid_t(xy) ),
    grid_( unstructured_grid(get()) ) {
}


inline const StructuredGrid::grid_t* structured_grid( const Grid::Implementation *grid ) {
    return dynamic_cast<const StructuredGrid::grid_t*>(grid);
}

StructuredGrid::StructuredGrid():
    Grid(),
    grid_( nullptr ) {
}

StructuredGrid::StructuredGrid( const Grid& grid ):
    Grid( grid ),
    grid_( structured_grid(get()) ) {
}

StructuredGrid::StructuredGrid( const Grid::Implementation* grid ):
    Grid( grid ),
    grid_( structured_grid(get()) ) {
}

StructuredGrid::StructuredGrid( const std::string& grid, const Domain& domain ):
    Grid( grid, domain ),
    grid_( structured_grid(get()) ) {
}

StructuredGrid::StructuredGrid( const Config& grid ):
    Grid( grid ),
    grid_( structured_grid(get()) ) {
}

StructuredGrid::StructuredGrid(
    const XSpace&     xspace,
    const YSpace&     yspace,
    const Projection& projection,
    const Domain&     domain
) :
    Grid( new detail::grid::Structured( xspace, yspace, projection, domain ) ),
    grid_( structured_grid(get()) ) {
}


ReducedGaussianGrid::ReducedGaussianGrid( const std::vector<long>& nx, const Domain& domain ):
    ReducedGaussianGrid::grid_t( detail::grid::reduced_gaussian(nx,domain) ) {
}

ReducedGaussianGrid::ReducedGaussianGrid( const std::initializer_list<long>& nx ):
    ReducedGaussianGrid( std::vector<long>(nx) ) {
}


} // namespace grid
} // namespace atlas
