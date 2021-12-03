/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <initializer_list>

#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/grid/Unstructured.h"


namespace atlas {

//---------------------------------------------------------------------------------------------------------------------
// Further grid interpretation classes defined in this file

class UnstructuredGrid;

/*
                                             Grid
                                               |
                                    +----------+----------+
                                    |                     |
                             StructuredGrid        UnstructuredGrid

*/

//---------------------------------------------------------------------------------------------------------------------

class UnstructuredGrid : public Grid {
public:
    using grid_t = grid::detail::grid::Unstructured;

public:
    UnstructuredGrid();
    UnstructuredGrid(const Grid&);
    UnstructuredGrid(const Config&);
    UnstructuredGrid(const Grid::Implementation*);
    UnstructuredGrid(std::vector<PointXY>*);  // takes ownership
    UnstructuredGrid(std::initializer_list<PointXY>);
    UnstructuredGrid(const Grid&, const Domain&);  // Create a new unstructured grid!

    operator bool() const { return valid(); }
    UnstructuredGrid(std::vector<PointXY>&&);       // move constructor
    UnstructuredGrid(const std::vector<PointXY>&);  // creates copy

    bool valid() const { return grid_; }

    using Grid::lonlat;
    using Grid::xy;

    // @brief Copy coordinates x,y of point `n` into xy[2]
    void xy(idx_t n, double xy[]) const { grid_->xy(n, xy); }

    // @brief Copy coordinates lon,lat of point `n` into lonlat[2]
    void lonlat(idx_t n, double lonlat[]) const { grid_->lonlat(n, lonlat); }

    // @brief Create PointXY at grid index `n`
    PointXY xy(idx_t n) const { return grid_->xy(n); }

    // @brief Create PointLonLat at grid index `n`
    PointLonLat lonlat(idx_t n) const { return grid_->lonlat(n); }

private:
    const grid_t* grid_;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
