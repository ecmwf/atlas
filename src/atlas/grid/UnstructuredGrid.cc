/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/UnstructuredGrid.h"

#include <limits>
#include <string>
#include <vector>

#include "eckit/config/Parametrisation.h"

#include "atlas/domain/Domain.h"
#include "atlas/grid/Spacing.h"
#include "atlas/projection/Projection.h"
#include "atlas/util/Config.h"

namespace atlas {

inline const UnstructuredGrid::grid_t* unstructured_grid( const Grid::Implementation* grid ) {
    return dynamic_cast<const UnstructuredGrid::grid_t*>( grid );
}

UnstructuredGrid::UnstructuredGrid() : Grid() {}

UnstructuredGrid::UnstructuredGrid( const Grid& grid ) : Grid( grid ), grid_( unstructured_grid( get() ) ) {}

UnstructuredGrid::UnstructuredGrid( const Grid::Implementation* grid ) :
    Grid( grid ), grid_( unstructured_grid( get() ) ) {}

UnstructuredGrid::UnstructuredGrid( const Config& grid ) : Grid( grid ), grid_( unstructured_grid( get() ) ) {}

UnstructuredGrid::UnstructuredGrid( std::vector<PointXY>* xy ) :
    Grid( new UnstructuredGrid::grid_t( xy ) ), grid_( unstructured_grid( get() ) ) {}

UnstructuredGrid::UnstructuredGrid( std::vector<PointXY>&& xy ) :
    Grid( new UnstructuredGrid::grid_t( std::forward<std::vector<PointXY>>( xy ) ) ),
    grid_( unstructured_grid( get() ) ) {}

UnstructuredGrid::UnstructuredGrid( const std::vector<PointXY>& xy ) :
    Grid( new UnstructuredGrid::grid_t( std::forward<const std::vector<PointXY>>( xy ) ) ),
    grid_( unstructured_grid( get() ) ) {}

UnstructuredGrid::UnstructuredGrid( std::initializer_list<PointXY> xy ) :
    Grid( new UnstructuredGrid::grid_t( xy ) ), grid_( unstructured_grid( get() ) ) {}

UnstructuredGrid::UnstructuredGrid( const Grid& grid, const Grid::Domain& domain ) :
    Grid( new UnstructuredGrid::grid_t( *grid.get(), domain ) ), grid_( unstructured_grid( get() ) ) {}


}  // namespace atlas
