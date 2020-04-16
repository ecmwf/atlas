/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <functional>
#include <initializer_list>
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
                                    +----------+----------+------------------------+
                                    |                     |                        |
                             StructuredGrid        UnstructuredGrid           CubedSphere
                                    |                                              |
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

    inline idx_t ny() const { return grid_->ny(); }

    inline idx_t nx( idx_t j ) const { return grid_->nx( j ); }

    inline const std::vector<idx_t>& nx() const { return grid_->nx(); }

    //inline const std::vector<double>& y() const { return grid_->y(); }

    /// x coordinate for given grid point {i,j}
    inline double x( idx_t i, idx_t j, idx_t t ) const { return grid_->x( i, j, t ); }

    /// y coordinate for given grid row {j}
    inline double y( idx_t i, idx_t j, idx_t t ) const { return grid_->y( i, j, t ); }

    using Grid::xy;
    void xy( idx_t i, idx_t j, idx_t t, double xy[] ) const { grid_->xy( i, j, t, xy ); }

    using Grid::lonlat;
    void lonlat( idx_t i, idx_t j, idx_t t, double lonlat[] ) const { grid_->lonlat( i, j, t, lonlat ); }

    PointXY xy( idx_t i, idx_t j, idx_t t) const { return PointXY( x( i, j, t), y( i, j, t ) ); }

    PointLonLat lonlat( idx_t i, idx_t j, idx_t t ) const { return grid_->lonlat( i, j, t ); }

    inline int GetCubeNx() const { return grid_->GetCubeNx(); }

private:
    const grid_t* grid_;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
