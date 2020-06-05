/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/grid/CubedSphere.h"

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------
class CubedSphereGrid;
class EquiDistCubedSphereGrid;
class EquiAnglCubedSphereGrid;

/*
                                             Grid
                                               |
                                    +----------+----------+-------------------------------+
                                    |                     |                               |
                             StructuredGrid        UnstructuredGrid                  CubedSphere
                                    |                                                     |
               +--------------------+-----------------------+               +------+------+------+------+
               |                    |                       |               |             |             |
          ReducedGrid          GaussianGrid            RegularGrid       EquiDist      EquiAngl     EquiDistFV3
               |                 |     |                 |     |
               +--------+--------+     +--------+--------+     +-----+
                        |                       |                    |
               ReducedGaussianGrid     RegularGaussianGrid    RegularLonLatGrid
*/

//---------------------------------------------------------------------------------------------------------------------

/*! @class CubedSphereGrid
 * @brief Specialization of Grid, where the grid is a cubed sphere
 * @details
 * @copydetails Grid
 * @dotfile classatlas_1_1CubedSphereGrid__inherit__graph.dot
 */


/*
Description of indices
----------------------
This class and the implementations in detail/grid/CubedSphere.h are for the structured cubed-sphere
grid.

In the following there are indices in the grid, denoted idx_t i, idx_t j, idx_t t and positions in
the grid, denoted double x, double y and t. For the cubed-sphere each face of the tile is
represented by rank 2 arrays with index in that array denoted by i and j. The face is denoted by t.

For a regular like grid, translating between i, j and x, y could be done from a starting point and
the grid spacing. But this is not possible for the cubed-sphere grid, where the spacing is not
necessarily regular and the array not rank 2 for the entire grid, i.e. is made more complicated
by the tile index.

x and y represent the position on a plane with the tiles arranged in an 'unfolded manner', shown
below. Tile t is carried around also as part of the array xy in order to do the projection properly.

Cubed-sphere Panel Arrangement:

y

^             .......
|            |       :
|            |   3   :
|            |       :
|             --->---
|   *.......  .......  -------  -------
|   |       :|       :|       :|       :
|   |   1   :|   2   :v   4   :v   5   :
|   |       :|       :|       :|       :
|     --->---  --->--* .......  .......
|                               -------
|                              |       :
|                              v   6   :
|                              |       :
|                               .......
----------------------------------------------------------->  x

  Key
  Solid lines: left and bottom edges of panel
  Dotted lines: right and top edges of panel
  >, <, v: direction of increasing index in first dimension
  *: location of two extra points (ngrid = 6 * (NCube+1) * (NCube+1) + 2)

*/

class CubedSphereGrid : public Grid {
public:
    using grid_t = grid::detail::grid::CubedSphere;

public:
  CubedSphereGrid();
  CubedSphereGrid( const Grid& );
  CubedSphereGrid( const Grid::Implementation* );
  CubedSphereGrid( const std::string& name );
  CubedSphereGrid( const Config& );
  CubedSphereGrid( const int&, const Projection& = Projection() );

  operator bool() const { return valid(); }

  bool valid() const { return grid_; }

  using Grid::xy;
  // Given indexes in the array (i, j, t) return position array xyt
  void xy( idx_t i, idx_t j, idx_t t, double xyt[] ) const { grid_->xy( i, j, t, xyt ); }

  // Given indexes in the array (i, j, t) return position array as a PointXY object.
  PointXY xy( idx_t i, idx_t j, idx_t t) const { return grid_->xy( i, j, t ); }

  using Grid::lonlat;
  // Given indexes in the array (i, j, t) return lat/lon array (via the projection)
  void lonlat( idx_t i, idx_t j, idx_t t, double lonlat[] ) const { grid_->lonlat( i, j, t, lonlat ); }

  // Given indexes in the array (i, j, t) return lat/lon as a PointLonLat object
  PointLonLat lonlat( idx_t i, idx_t j, idx_t t ) const { return grid_->lonlat( i, j, t ); }

  // Convect lonlat array into indices array (via inverse projection)
  void lonlat2xy( double lonlat[], idx_t ijt[] ) const { grid_->lonlat2xy( lonlat, ijt ); }

  // Return the size of the cubed sphere grid, where CubeNX is the number of grid boxes along the edge of a tile
  inline const int GetCubeNx() const { return grid_->GetCubeNx(); }

  // Return the number of tiles
  inline const int GetNTiles() const { return grid_->GetNTiles(); }

private:
  const grid_t* grid_;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
