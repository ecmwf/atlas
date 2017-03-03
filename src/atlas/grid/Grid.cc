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
#include "atlas/grid/Projection.h"
#include "atlas/util/Config.h"
#include "atlas/grid/detail/grid/Gaussian.h"

namespace atlas {
namespace grid {

Grid::Grid():
    grid_( nullptr ) {
}

Grid::Grid(const Grid& grid):
    grid_( grid.grid_ ) {
}

Grid::Grid( const Grid::grid_t *grid ):
    grid_( grid ) {
}

Grid::Grid( const std::string& shortname ) {
    grid_ = Grid::grid_t::create( shortname );
}

Grid::Grid( const Config& p ) {
    grid_ = Grid::grid_t::create(p);
}

const StructuredGrid::grid_t* structured_grid( const Grid::grid_t *grid ) {
    const StructuredGrid::grid_t* g( dynamic_cast<const StructuredGrid::grid_t*>(grid) );
    return g;
}

StructuredGrid::StructuredGrid():
    Grid(),
    grid_( nullptr ) {
}

StructuredGrid::StructuredGrid( const Grid& grid ):
    Grid( grid ),
    grid_( structured_grid(get()) ) {
}

StructuredGrid::StructuredGrid( const Grid::grid_t* grid ):
    Grid( grid ),
    grid_( structured_grid(get()) ) {
}

StructuredGrid::StructuredGrid( const std::string& grid ):
    Grid( grid ),
    grid_( structured_grid(get()) ) {
}

StructuredGrid::StructuredGrid( const Config& grid ):
    Grid( grid ),
    grid_( structured_grid(get()) ) {
}

const RegularGrid::grid_t* regular_grid( const Grid::grid_t* grid ) {
    const RegularGrid::grid_t* g( dynamic_cast<const RegularGrid::grid_t*>(grid) );
    if( g && g->reduced() ) {
        return nullptr;
    }
    return g;
}

RegularGrid::grid_t* RegularGrid::create( const Config& ) {
  NOTIMP;
  // return nullptr;
}

RegularGrid::RegularGrid():
    StructuredGrid(),
    grid_( nullptr ),
    nx_(0) {
}

RegularGrid::RegularGrid( const Grid& grid ):
    StructuredGrid( grid ),
    grid_( regular_grid(get()) ) {
    if( grid_ ) nx_ = StructuredGrid::nx().front();
}

RegularGrid::RegularGrid( const detail::grid::Grid *grid ):
    StructuredGrid(grid),
    grid_( regular_grid(get()) ) {
    if( grid_ ) nx_ = StructuredGrid::nx().front();
}

RegularGrid::RegularGrid( const std::string& grid ):
    StructuredGrid(grid),
    grid_( regular_grid(get()) ) {
    if( grid_ ) nx_ = StructuredGrid::nx().front();
}


RegularGrid::RegularGrid( const Config& p ):
    StructuredGrid(create(p)),
    grid_( regular_grid(get()) ) {
    if( grid_ ) nx_ = StructuredGrid::nx().front();
}

ReducedGaussianGrid::ReducedGaussianGrid( const std::vector<long>& nx ):
    ReducedGaussianGrid::Grid( detail::grid::reduced_gaussian(nx) ) {
}

ReducedGaussianGrid::ReducedGaussianGrid( const std::initializer_list<long>& nx ):
    ReducedGaussianGrid( std::vector<long>(nx) ) {
}


} // namespace Grid
} // namespace atlas
