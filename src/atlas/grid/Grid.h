/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include <string>

#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/grid/detail/grid/Regular.h"
#include "atlas/grid/Projection.h"
#include "atlas/grid/Domain.h"
#include "atlas/grid/Iterator.h"

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------
// classes defined in this file

class Grid;
class StructuredGrid;
class RegularGrid;
class GaussianGrid;
class RegularGaussianGrid;
class ReducedGaussianGrid;
class RegularLonLatGrid;
class ShiftedLonLatGrid;

//---------------------------------------------------------------------------------------------------------------------

class Grid {

public:

    using grid_t     = detail::grid::Grid;
    using Config     = grid_t::Config;
    using Spec       = grid_t::Spec;
    using Domain     = grid::Domain;
    using Projection = grid::Projection;
    using Iterator   = grid::Iterator;
    
    using iterator = Iterator;
    using const_iterator = iterator;

public:

    Grid();
    Grid( const Grid& );
    Grid( const grid_t* );
    Grid( const std::string& );
    Grid( const Config& );

    operator bool() const { return grid_; }
    operator const grid_t&() const { return *get(); }
    bool operator==( const Grid& other ) const { return grid_ == other.grid_; }
    bool operator!=( const Grid& other ) const { return grid_ != other.grid_; }

    size_t npts() const { return grid_->npts(); }

    const Projection& projection() const { return grid_->projection(); }
    const Domain& domain() const { return grid_->domain(); }
    std::string name() const { return grid_->name(); }
    std::string uid() const { return grid_->uid(); }

    Spec spec() const { return grid_->spec(); }
    
    iterator begin() const { return grid_->begin(); }
    iterator end()   const { return grid_->end(); }

    const grid_t* get() const { return grid_.get(); }

private:

    eckit::SharedPtr<const grid_t> grid_;

};

//---------------------------------------------------------------------------------------------------------------------

class StructuredGrid: public Grid {

public:

  using grid_t = detail::grid::Structured;
  using YSpace = grid_t::YSpace;

public:

    StructuredGrid();
    StructuredGrid( const Grid& );
    StructuredGrid( const Grid::grid_t* );
    StructuredGrid( const std::string& );
    StructuredGrid( const Config& );

    operator bool() const {
      return grid_ && valid();
    }

    inline size_t ny() const {
        return grid_->ny();
    }

    inline size_t nx( size_t j ) const {
        return grid_->nlon(j);;
    }

    inline const std::vector<long>& nx() const {
        return grid_->pl();
    }

    inline const long nxmax() const {
        return grid_->nlonmax();
    }

    inline const std::vector<double>& y() const {
        return grid_->latitudes();
    }

    inline double x( const size_t i, const size_t j ) const {
        return grid_->lon(j,i);
    }

    inline double y( const size_t j ) const {
        return grid_->lat(j);
    }

    void xy( const size_t i, const size_t j, double xy[] ) const {
      xy[0] = x(i,j);
      xy[1] = y(j) ;
    }

    PointXY xy( const size_t i, const size_t j ) const {
      return PointXY( x(i,j), y(j) );
    }

    PointLonLat lonlat( const size_t i, const size_t j ) const {
      return grid_->geolonlat(i,j);
    }

    inline bool reduced() const {
        return grid_->reduced();
    }

    bool periodic() const {
        return grid_->isPeriodicX();
    }

    const YSpace& yspace() const {
      return grid_->yspace();
    }

private:

    virtual bool valid() const { return true; }
    const grid_t* grid_;
};

//---------------------------------------------------------------------------------------------------------------------

class RegularGrid: public StructuredGrid {

public:

    using grid_t = StructuredGrid::grid_t;

public:

    RegularGrid();
    RegularGrid( const Grid& );
    RegularGrid( const Grid::grid_t* );
    RegularGrid( const std::string& );
    RegularGrid( const Config& );

    operator bool() const { return grid_ && valid(); }

    size_t nx() const { return nx_; }

    inline double x( const size_t i ) const {
        return StructuredGrid::x(i,0);
    }

    using StructuredGrid::y;

    inline double y( const size_t j ) const {
        return grid_->lat(j);
    }

    using StructuredGrid::xy;

    PointXY xy( const size_t i, const size_t j ) const {
      return PointXY( x(i), y(j ) );
    }

private:

    virtual bool valid() const override { return true; }
    static grid_t* create( const Config& );
    const grid_t* grid_;
    size_t nx_ = {0};
};

//---------------------------------------------------------------------------------------------------------------------

template< class Grid >
class Gaussian : public Grid {

public:

    using Grid::Grid;

    long N() const { return Grid::ny()/2; }

    inline double lon( const size_t i, const size_t j ) const {
        return Grid::x(i,j);
    }

    inline double lat( const size_t j ) const {
        return Grid::y(j);
    }

    PointLonLat lonlat( const size_t i, const size_t j ) const {
      return Grid::xy(i,j);
    }

protected:

    bool gaussian() const {
      return Grid::domain().global()
        &&   not Grid::projection()
        &&   Grid::yspace().type() == "gaussian";
    }
};

//---------------------------------------------------------------------------------------------------------------------

class LonLat : public RegularGrid {

public:

    using RegularGrid::RegularGrid;

public:

  inline double lon( const size_t i ) const {
      return x(i);
  }

  inline double lat( const size_t j ) const {
      return y(j);
  }

  PointLonLat lonlat( const size_t i, const size_t j ) const {
    return xy(i,j);
  }

protected:

    bool global_lonlat() const {
      return domain().global()
        &&   not projection()
        &&   yspace().type() == "linear";
    }

    bool regular_lon() const {
        return x(0) == 0.;
    }

    bool regular_lat() const {
        return y(0) == 90.
            && ny()%2 == 1;
    }

    bool shifted_lon() const {
       return   x(0) == 0.5*360./nx();
    }

    bool shifted_lat() const {
        return y(0) == 90.-0.5*180./ny()
            && ny()%2 == 0;
    }

private:

    virtual bool valid() const =0;
};

//---------------------------------------------------------------------------------------------------------------------

class GaussianGrid : public Gaussian<StructuredGrid> {

    using Grid = Gaussian<StructuredGrid>;

public:

    using Grid::Grid;

private:

    virtual bool valid() const override {
        return gaussian();
    }
};

//---------------------------------------------------------------------------------------------------------------------

class ReducedGaussianGrid : public Gaussian<StructuredGrid> {

    using Grid = Gaussian<StructuredGrid>;

public:

    using Grid::Grid;

private:

    virtual bool valid() const override {
        return gaussian() && reduced();
    }
};

//---------------------------------------------------------------------------------------------------------------------

class RegularGaussianGrid : public Gaussian<RegularGrid> {

    using Grid = Gaussian<RegularGrid>;

public:

    using Grid::Grid;

    inline double lon( const size_t i ) const {
        return x(i);
    }

    inline double lat( const size_t j ) const {
        return y(j);
    }

    PointLonLat lonlat( const size_t i, const size_t j ) const {
      return xy(i,j);
    }

private:

    virtual bool valid() const override {
        return gaussian();
    }
};

//---------------------------------------------------------------------------------------------------------------------

class RegularLonLatGrid : public LonLat {

public:

    using LonLat::LonLat;

private:

    virtual bool valid() const override {
        return global_lonlat()
            && regular_lon()
            && regular_lat();
    }
};

//---------------------------------------------------------------------------------------------------------------------

class ShiftedLonLatGrid : public LonLat {

public:

    using LonLat::LonLat;

private:

    virtual bool valid() const override {
        return global_lonlat()
            && shifted_lon()
            && shifted_lat();
    }
};

//---------------------------------------------------------------------------------------------------------------------

}
}
