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
               +--------------------+-----------------------+               +------+------+
               |                    |                       |               |             |
          ReducedGrid          GaussianGrid            RegularGrid       EquiDist      EquiAngl
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
|            |   2   :
|            |       :
|             --->---
|   *.......  .......  -------  -------
|   |       :|       :|       :|       :
|   |   0   :|   1   :v   3   :v   4   :
|   |       :|       :|       :|       :
|     --->---  --->--* .......  .......
|                               -------
|                              |       :
|                              v   5   :
|                              |       :
|                               .......
----------------------------------------------------------->  x

Key
Solid lines: left and bottom edges of panel
Dotted lines: right and top edges of panel
>, <, v: direction of increasing index in first dimension
*: location of two extra points (ngrid = 6 * NCube * NCube + 2)

----------------------------------------------------------------------------------------------------

NCube = 3 example ( 6 * NCube * NCube + 2 = 56 grid points )

Index of each point within the grid (i,j,t):

  () denotes actual points in the grids
  <> denotes points that are duplicates (not part of the grid) but shown for completeness of the grid



                        <0,3,0>-<0,2,4>-<0,1,4>-<0,0,4>
                           |                       |
                           |                       |
                        (0,1,2) (1,2,2) (2,2,2) <0,2,3>
                           |                       |
                           |                       |
                        (0,1,2) (1,1,2) (2,1,2) <0,1,3>
                           |                       |
                           |                       |
(0,3,0)-<0,1,2>-<0,1,2>-(0,0,2)-(1,0,2)-(2,0,2)-(0,0,3)-(0,1,3)-(0,2,3)-(0,0,4)-(0,1,4)-(0,2,4)-<0,3,0>
  |                        |                       |                       |                       |
  |                        |                       |                       |                       |
(0,2,0) (1,2,0) (2,2,0) (0,1,1) (1,2,1) (2,2,1) (1,0,3) (1,1,3) (1,2,3) (1,0,4) (1,1,4) (1,2,4) <0,1,0>
  |                        |                       |                       |                       |
  |                        |                       |                       |                       |
(0,1,0) (1,1,0) (2,1,0) (0,1,1) (1,1,1) (2,1,1) (2,0,3) (2,1,3) (2,2,3) (2,0,4) (2,1,4) (2,2,4) <0,1,0>
  |                        |                       |                       |                       |
  |                        |                       |                       |                       |
(0,0,0)-(1,0,0)-(2,0,0)-(0,0,1)-(1,0,1)-(2,0,1)-(3,0,1)-<2,0,5>-<1,0,5>-(0,0,5)-(0,1,5)-(0,2,5)-<0,0,0>
                                                                           |                       |
                                                                           |                       |
                                                                        (1,0,5) (1,1,5) (1,2,5) <1,0,0>
                                                                           |                       |
                                                                           |                       |
                                                                        (2,0,5) (2,1,5) (2,2,5) <2,0,0>
                                                                           |                       |
                                                                           |                       |
                                                                        <3,0,1>-<2,0,1>-<1,0,1>-<0,0,1>



Position of each point within the grid (x,y,t):

  () denotes actual points in the grids
  +  denotes duplicate points

  For the xyt, t is only required in order to pass to the projection and thus know how to rotate
  the projection of tile 0 to obtain the other tiles


                           +-------+-------+-------+
                           |                       |
                           |                       |
                        (3,8,2) (4,8,2) (5,8,2)    +
                           |                       |
                           |                       |
                        (3,7,2) (4,7,2) (5,7,2)    +
                           |                       |
                           |                       |
(0,6,0)----+-------+----(3,6,2)-(4,6,2)-(5,6,2)-(6,6,3)-(7,6,3)-(8,6,3)-(9,6,4)-(10,6,4)-(11,6,4)----+
   |                       |                       |                       |                         |
   |                       |                       |                       |                         |
(0,5,0) (1,5,0) (2,5,0) (3,5,1) (4,5,1) (5,5,1) (6,5,3) (7,5,3) (8,5,3) (9,5,4) (10,5,4) (11,5,4)    +
   |                       |                       |                       |                         |
   |                       |                       |                       |                         |
(0,4,0) (1,4,0) (2,4,0) (3,4,1) (4,4,1) (5,4,1) (6,4,3) (7,4,3) (8,4,3) (9,4,4) (10,4,4) (11,4,4)    +
   |                       |                       |                       |                         |
   |                       |                       |                       |                         |
(0,3,0)-(1,3,0)-(2,3,0)-(3,3,1)-(4,3,1)-(5,3,1)-(6,3,1)----+-------+----(9,3,5)-(10,3,5)-(11,3,5)----+
                                                                           |                         |
                                                                           |                         |
                                                                        (9,2,5) (10,2,5) (11,2,5)    +
                                                                           |                         |
                                                                           |                         |
                                                                        (9,1,5) (10,1,5) (11,1,5)    +
                                                                           |                         |
                                                                           |                         |
                                                                           +--------+-------+--------+


Position in the grid iterator/mesh:

  First element of each tile: {0, 10, 20, 29, 38, 47}

             +---+---+---+
             |           |
            26  27  28   +
             |           |
            23  24  25   +
             |           |
 9---+---+--20--21--22--29--32--35---38--41--44---+
 |           |           |            |           |
 6   7   8  17  18  19  30  33  36   39  42  45   +
 |           |           |            |           |
 3   4   5  14  15  16  31  34  37   40  43  46   +
 |           |           |            |           |
 0---1---2--10--11--12--13---+---+---47--50--53---+
                                      |           |
                                     48  51  54   +
                                      |           |
                                     49  52  55   +
                                      |           |
                                      +---+---+---+

----------------------------------------------------------------------------------------------------

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
  void xyt( idx_t i, idx_t j, idx_t t, double xyt[] ) const { grid_->xyt( i, j, t, xyt ); }
  PointXY xyt( idx_t i, idx_t j, idx_t t) const { return grid_->xyt( i, j, t ); }
  // Given indexes in the array (i, j, t) return position array xyt

  void xy( idx_t i, idx_t j, idx_t t, double xy[] ) const { grid_->xy( i, j, t, xy ); }
  PointXY xy( idx_t i, idx_t j, idx_t t) const { return grid_->xy( i, j, t); }

  using Grid::lonlat;
  // Given indexes in the array (i, j) return lat/lon array (via the projection)
  void lonlat( idx_t i, idx_t j, idx_t t, double lonlat[] ) const { grid_->lonlat( i, j, t, lonlat ); }

  // Given indexes in the array (i, j, t) return lat/lon as a PointLonLat object
  PointLonLat lonlat( idx_t i, idx_t j, idx_t t ) const { return grid_->lonlat( i, j, t ); }

  // Return the size of the cubed sphere grid, where CubeNX is the number of grid boxes along the edge of a tile
  inline int GetCubeNx() const { return grid_->GetCubeNx(); }

  // Return the number of tiles
  inline int GetNTiles() const { return grid_->GetNTiles(); }

  // Transform from xy space to xyt space that is a function of resolution.
  void xy2xyt(const double xy[], double xyt[]) const {return grid_->xy2xyt(xy, xyt); }

  // Transform from xyt space to xy space
  void xyt2xy(const double xyt[], double xy[]) const {return grid_->xyt2xy(xyt, xy); }

private:
  const grid_t* grid_;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas