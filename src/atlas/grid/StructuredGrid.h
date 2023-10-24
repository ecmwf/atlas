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

#include <functional>
#include <initializer_list>
#include <string>

#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/grid/Healpix.h"
#include "atlas/grid/detail/grid/Structured.h"

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------
class StructuredGrid;
class RegularGrid;
class GaussianGrid;
class ReducedGaussianGrid;
class RegularGaussianGrid;
class RegularLonLatGrid;
class ShiftedLonLatGrid;
class HealpixGrid;

/*
                                                 Grid
												   |
                                        +----------+----------+
                                        |                     |
                                  StructuredGrid       UnstructuredGrid
                                        |
                   +--------------------+-----------------------+------------------+
                   |                    |                       |                  |
              ReducedGrid          GaussianGrid            RegularGrid        HealpixGrid
                   |                 |     |                 |     |
                   +--------+--------+     +--------+--------+     +-----+
                            |                       |                    |
                   ReducedGaussianGrid     RegularGaussianGrid    RegularLonLatGrid
*/

//---------------------------------------------------------------------------------------------------------------------

/*! @class StructuredGrid
 * @brief Specialization of Grid, where the grid can be represented by rows with uniform distribution
 * @details 
 * @copydetails Grid
 * @dotfile classatlas_1_1StructuredGrid__inherit__graph.dot
 */
class StructuredGrid : public Grid {
public:
    using grid_t = grid::detail::grid::Structured;
    using XSpace = grid_t::XSpace;
    using YSpace = grid_t::YSpace;

public:
    StructuredGrid();
    StructuredGrid(const Grid&);
    StructuredGrid(const Grid::Implementation*);
    StructuredGrid(const std::string& name, const Domain& = Domain());
    StructuredGrid(const std::string& name, const Projection&, const Domain& = Domain());
    StructuredGrid(const Config&);
    StructuredGrid(const XSpace&, const YSpace&, const Projection& = Projection(), const Domain& = Domain());
    StructuredGrid(const Grid&, const Domain&);

    operator bool() const { return valid(); }

    bool valid() const { return grid_; }

    inline idx_t ny() const { return grid_->ny(); }

    inline idx_t nx(idx_t j) const { return grid_->nx(j); }

    inline const std::vector<idx_t>& nx() const { return grid_->nx(); }

    inline idx_t nxmax() const { return grid_->nxmax(); }

    inline const std::vector<double>& y() const { return grid_->y(); }

    /// x coordinate for given grid point {i,j}
    inline double x(idx_t i, idx_t j) const { return grid_->x(i, j); }

    /// y coordinate for given grid row {j}
    inline double y(idx_t j) const { return grid_->y(j); }

    /// increment in x for a given grid row {j}
    inline double dx(idx_t j) const { return grid_->dx(j); }

    /// x coordinate of beginning of a given grid row {j}
    inline double xmin(idx_t j) const { return grid_->xmin(j); }


    using Grid::xy;
    void xy(idx_t i, idx_t j, double xy[]) const { grid_->xy(i, j, xy); }

    using Grid::lonlat;
    void lonlat(idx_t i, idx_t j, double lonlat[]) const { grid_->lonlat(i, j, lonlat); }

    PointXY xy(idx_t i, idx_t j) const { return PointXY(x(i, j), y(j)); }

    PointLonLat lonlat(idx_t i, idx_t j) const { return grid_->lonlat(i, j); }

    inline bool reduced() const { return grid_->reduced(); }

    inline bool regular() const { return not reduced(); }

    bool periodic() const { return grid_->periodic(); }

    const XSpace& xspace() const { return grid_->xspace(); }

    const YSpace& yspace() const { return grid_->yspace(); }

    gidx_t index(idx_t i, idx_t j) const { return grid_->index(i, j); }

    void index2ij(gidx_t gidx, idx_t& i, idx_t& j) const { grid_->index2ij(gidx, i, j); }

private:
    const grid_t* grid_;
};

//---------------------------------------------------------------------------------------------------------------------

/// @class ReducedGrid
/// @brief Specialization of StructuredGrid, where not all rows have the same number of grid points
class ReducedGrid : public StructuredGrid {
public:
    using StructuredGrid::StructuredGrid;

    operator bool() const { return valid(); }

    bool valid() const { return StructuredGrid::valid() && reduced(); }
};

//---------------------------------------------------------------------------------------------------------------------

/// @class RegularGrid
/// @brief Specialization of StructuredGrid, where all rows have the same number of grid points
class RegularGrid : public StructuredGrid {
public:
    using StructuredGrid::dx;
    using StructuredGrid::StructuredGrid;
    using StructuredGrid::x;
    using StructuredGrid::xy;

    operator bool() const { return valid(); }

    bool valid() const { return StructuredGrid::valid() && regular(); }

    idx_t nx() const { return nxmax(); }

    inline double dx() const { return dx(0); }

    inline double x(idx_t i) const { return x(i, 0); }

    PointXY xy(idx_t i, idx_t j) const { return PointXY(x(i), y(j)); }
};

//---------------------------------------------------------------------------------------------------------------------

template <class Grid>
class Gaussian : public Grid {
public:
    using Grid::Grid;

    idx_t N() const { return Grid::ny() / 2; }

protected:
    bool gaussian() const { return Grid::domain().global() && Grid::yspace().type() == "gaussian"; }
};

//---------------------------------------------------------------------------------------------------------------------

/// @class GaussianGrid
/// @brief Specialization of StructuredGrid, where rows follow a Gaussian distribution
class GaussianGrid : public Gaussian<StructuredGrid> {
    using grid_t = Gaussian<StructuredGrid>;

public:
    using grid_t::grid_t;

    operator bool() const { return valid(); }

    bool valid() const { return StructuredGrid::valid() && gaussian(); }
};

//---------------------------------------------------------------------------------------------------------------------

/// @class ReducedGaussianGrid
/// @brief Specialization of ReducedGrid, where rows follow a Gaussian distribution
class ReducedGaussianGrid : public Gaussian<ReducedGrid> {
    using grid_t = Gaussian<ReducedGrid>;

public:
    using grid_t::grid_t;
    ReducedGaussianGrid(const std::initializer_list<idx_t>& pl);
    ReducedGaussianGrid(const std::vector<long>& pl, const Domain& = Domain());
    ReducedGaussianGrid(const std::vector<int>& pl, const Domain& = Domain());
    ReducedGaussianGrid(const std::vector<long>& pl, const Projection&);
    ReducedGaussianGrid(const std::vector<int>& pl, const Projection&);

    operator bool() const { return valid(); }

    bool valid() const { return ReducedGrid::valid() && gaussian(); }
};

//---------------------------------------------------------------------------------------------------------------------

/// @class RegularGaussianGrid
/// @brief Specialization of RegularGaussianGrid, where rows follow a Gaussian distribution
class RegularGaussianGrid : public Gaussian<RegularGrid> {
    using grid_t = Gaussian<RegularGrid>;

public:
    using grid_t::grid_t;
    RegularGaussianGrid(int N, const Domain& = Domain());

    operator bool() const { return valid(); }

    bool valid() const { return RegularGrid::valid() && gaussian(); }
};

//---------------------------------------------------------------------------------------------------------------------

/// @class RegularLonLatGrid
/// @brief Specialization of RegularGrid, assuming a global domain
class RegularLonLatGrid : public RegularGrid {
public:
    using RegularGrid::RegularGrid;

public:
    operator bool() const { return valid(); }

    bool valid() const { return RegularGrid::valid() && global_lonlat(); }

    inline double lon(idx_t i) const { return x(i); }

    inline double lat(idx_t j) const { return y(j); }

    PointLonLat lonlat(idx_t i, idx_t j) const { return xy(i, j); }

    bool standard() const { return standard_lon() && standard_lat(); }
    bool shifted() const { return shifted_lon() && shifted_lat(); }
    bool shiftedLon() const { return shifted_lon() && standard_lat(); }
    bool shiftedLat() const { return standard_lon() && shifted_lat(); }

protected:
    bool global_lonlat() const { return domain().global() && not projection() && yspace().type() == "linear"; }

    bool standard_lon() const { return x(0) == 0.; }

    bool standard_lat() const { return y(0) == 90. && ny() % 2 == 1; }

    bool shifted_lon() const { return x(0) == 0.5 * 360. / nx(); }

    bool shifted_lat() const { return y(0) == 90. - 0.5 * 180. / ny() && ny() % 2 == 0; }
};

//---------------------------------------------------------------------------------------------------------------------

/// @class HealpixGrid
/// @brief Specialization of StructuredGrid, assuming a global domain
class HealpixGrid : public StructuredGrid {
public:
    using grid_t = grid::detail::grid::Healpix;

public:
    HealpixGrid(const Grid&);
    HealpixGrid(int N, const std::string& ordering = "ring");

    operator bool() const { return valid(); }

    bool valid() const { return grid_; }

    long N() const { return grid_->nxmax() / 4; }

private:
    const grid_t* grid_;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
