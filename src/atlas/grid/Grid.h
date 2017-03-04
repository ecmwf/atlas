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

#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/grid/Projection.h"
#include "atlas/grid/Domain.h"
#include "atlas/grid/Iterator.h"

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------
// grid classes defined in this file

class Grid;               
class StructuredGrid;     
class RegularGrid;        
class GaussianGrid;       
class ReducedGaussianGrid;
class RegularGaussianGrid;
class RegularLonLatGrid;  
class ShiftedLonLatGrid;  

/*
                                  Grid
                                    |
                                    |
                             StructuredGrid
                                    |
               +--------------------+-----------------------+
               |                    |                       |
          ReducedGrid          GaussianGrid            RegularGrid
               |                 |     |                 |     |
               +--------+--------+     +--------+--------+     +-----+
                        |                       |                    |
               ReducedGaussianGrid     RegularGaussianGrid       LonLatGrid
                                                                     |
                                                         +------------------------+
                                                         |                        |
                                                  RegularLonLatGrid        ShiftedLonLatGrid
*/

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
        return valid();
    }
    
    bool valid() const {
        return grid_;
    }

    inline size_t ny() const {
        return grid_->ny();
    }

    inline size_t nx( size_t j ) const {
        return grid_->nx(j);
    }

    inline const std::vector<long>& nx() const {
        return grid_->nx();
    }

    inline const long nxmax() const {
        return grid_->nxmax();
    }

    inline const std::vector<double>& y() const {
        return grid_->y();
    }

    inline double x( const size_t i, const size_t j ) const {
        return grid_->x(i,j);
    }

    inline double y( const size_t j ) const {
        return grid_->y(j);
    }

    void xy( const size_t i, const size_t j, double xy[] ) const {
        xy[0] = x(i,j);
        xy[1] = y(j) ;
    }

    PointXY xy( const size_t i, const size_t j ) const {
        return PointXY( x(i,j), y(j) );
    }

    PointLonLat lonlat( const size_t i, const size_t j ) const {
        return grid_->lonlat(i,j);
    }

    inline bool reduced() const {
        return grid_->reduced();
    }

    inline bool regular() const {
        return not reduced();
    }

    bool periodic() const {
        return grid_->periodic();
    }

    const YSpace& yspace() const {
        return grid_->yspace();
    }

private:

    const grid_t* grid_;
};

//---------------------------------------------------------------------------------------------------------------------

class ReducedGrid: public StructuredGrid {

public:

    using StructuredGrid::StructuredGrid;
    
    operator bool() const {
        return valid();
    }
    
    bool valid() const {
        return StructuredGrid::valid() && reduced();
    }

};

//---------------------------------------------------------------------------------------------------------------------

class RegularGrid: public StructuredGrid {

public:

    using StructuredGrid::StructuredGrid;
    using StructuredGrid::x;
    using StructuredGrid::xy;
    
    operator bool() const {
        return valid();
    }
    
    bool valid() const {
        return StructuredGrid::valid() && regular();
    }

    size_t nx() const {
        return nxmax();
    }

    inline double x( const size_t i ) const {
        return x(i,0);
    }

    PointXY xy( const size_t i, const size_t j ) const {
        return PointXY( x(i), y(j ) );
    }
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

class GaussianGrid : public Gaussian<StructuredGrid> {

    using Grid = Gaussian<StructuredGrid>;

public:

    using Grid::Grid;

    operator bool() const {
        return valid();
    }
    
    bool valid() const {
        return StructuredGrid::valid() && gaussian();
    }
};

//---------------------------------------------------------------------------------------------------------------------

class ReducedGaussianGrid : public Gaussian<ReducedGrid> {

    using Grid = Gaussian<ReducedGrid>;

public:

    using Grid::Grid;
    ReducedGaussianGrid( const std::initializer_list<long>& pl );
    ReducedGaussianGrid( const std::vector<long>& pl );

    operator bool() const {
        return valid();
    }
    
    bool valid() const {
        return ReducedGrid::valid() && gaussian();
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
    
    operator bool() const {
        return valid();
    }
    
    bool valid() const {
        return RegularGrid::valid() && gaussian();
    }

};

//---------------------------------------------------------------------------------------------------------------------

class LonLatGrid : public RegularGrid {

public:

    using RegularGrid::RegularGrid;

public:

    operator bool() const {
        return valid();
    }
  
    bool valid() const {
        return RegularGrid::valid() && global_lonlat();
    }

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
};

//---------------------------------------------------------------------------------------------------------------------

class RegularLonLatGrid : public LonLatGrid {

public:

    using LonLatGrid::LonLatGrid;

    operator bool() const {
        return valid();
    }
  
    bool valid() const {
        return LonLatGrid::valid()
            && regular_lon()
            && regular_lat();
    }
};

//---------------------------------------------------------------------------------------------------------------------

class ShiftedLonLatGrid : public LonLatGrid {

public:

    using LonLatGrid::LonLatGrid;

    operator bool() const {
        return valid();
    }
  
    bool valid() const {
        return LonLatGrid::valid()
            && shifted_lon()
            && shifted_lat();
    }
};

//---------------------------------------------------------------------------------------------------------------------

}
}
