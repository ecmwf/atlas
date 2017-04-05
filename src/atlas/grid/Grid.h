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

#include "eckit/memory/SharedPtr.h"
#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/grid/detail/grid/Unstructured.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/grid/Projection.h"
#include "atlas/grid/Domain.h"
#include "atlas/grid/Iterator.h"

namespace eckit {
  class MD5;
}

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------
// grid classes defined in this file

class Grid;
class UnstructuredGrid;
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
                                    +----------+----------+
                                    |                     |
                             StructuredGrid        UnstructuredGrid
                                    |
               +--------------------+-----------------------+
               |                    |                       |
          ReducedGrid          GaussianGrid            RegularGrid
               |                 |     |                 |     |
               +--------+--------+     +--------+--------+     +-----+
                        |                       |                    |
               ReducedGaussianGrid     RegularGaussianGrid    RegularLonLatGrid
*/

//---------------------------------------------------------------------------------------------------------------------

class Grid {

public:

    using grid_t     = detail::grid::Grid;
    using Config     = grid_t::Config;
    using Spec       = grid_t::Spec;
    using Domain     = grid::Domain;
    using Projection = grid::Projection;
    
    class IterateXY {
    public:
      using iterator       = grid::IteratorXY;
      using const_iterator = iterator;
    public:
      IterateXY(const grid_t& grid) : grid_(grid) {}
      iterator begin() const { return grid_.xy_begin(); }
      iterator end()   const { return grid_.xy_end(); }
    private:
      const grid_t& grid_;
    };
    
    class IterateLonLat {
    public:
      using iterator       = grid::IteratorLonLat;
      using const_iterator = iterator;
    public:
      IterateLonLat(const grid_t& grid) : grid_(grid) {}
      iterator begin() const { return grid_.lonlat_begin(); }
      iterator end()   const { return grid_.lonlat_end(); }
    private:
      const grid_t& grid_;
    };

public:
  
    IterateXY     xy()     const { return IterateXY(*grid_);     }
    IterateLonLat lonlat() const { return IterateLonLat(*grid_); }

    Grid();
    Grid( const Grid& );
    Grid( const grid_t* );
    Grid( const std::string& name, const Domain& = Domain() );
    Grid( const Config& );

    operator bool() const { return grid_; }
    operator const grid_t&() const { return *get(); }
    bool operator==( const Grid& other ) const { return grid_ == other.grid_; }
    bool operator!=( const Grid& other ) const { return grid_ != other.grid_; }

    size_t size() const { return grid_->size(); }

    const Projection& projection() const { return grid_->projection(); }
    const Domain& domain() const { return grid_->domain(); }
    std::string name() const { return grid_->name(); }
    std::string uid() const { return grid_->uid(); }
    
    /// Adds to the MD5 the information that makes this Grid unique
    void hash(eckit::MD5& md5) const { return grid_->hash(md5); }

    Spec spec() const { return grid_->spec(); }

    const grid_t* get() const { return grid_.get(); }

private:

    eckit::SharedPtr<const grid_t> grid_;

};

//---------------------------------------------------------------------------------------------------------------------

class GridLonLat {
  
  
};

//---------------------------------------------------------------------------------------------------------------------

class UnstructuredGrid: public Grid {

public:

  using grid_t = detail::grid::Unstructured;

public:

    UnstructuredGrid();
    UnstructuredGrid( const Grid& );
    UnstructuredGrid( const Config& );
    UnstructuredGrid( const Grid::grid_t* );
    UnstructuredGrid( std::vector<PointXY>* ); // takes ownership

    operator bool() const {
        return valid();
    }

    bool valid() const {
        return grid_;
    }

    using Grid::xy;
    void xy( const size_t n, double xy[] ) const {
      PointXY _xy = grid_->xy(n);
      xy[0] = _xy.x();
      xy[1] = _xy.y();
    }

    PointXY xy( const size_t n ) const {
        return grid_->xy(n);
    }

    PointLonLat lonlat( const size_t n ) const {
        return grid_->lonlat(n);
    }


private:

    const grid_t* grid_;
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
    StructuredGrid( const std::string& name, const Domain& = Domain() );
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

    using Grid::xy;
    void xy( const size_t i, const size_t j, double xy[] ) const {
        grid_->xy(i,j,xy);
    }

    void lonlat( const size_t i, const size_t j, double lonlat[] ) const {
        grid_->lonlat(i,j,lonlat);
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
    ReducedGaussianGrid( const std::vector<long>& pl, const Domain& domain = Domain() );

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

class RegularLonLatGrid : public RegularGrid {

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

    bool standard()      const { return standard_lon() && standard_lat(); }
    bool shifted()       const { return shifted_lon()  && shifted_lat();  }
    bool shiftedLon()    const { return shifted_lon()  && standard_lat(); }
    bool shiftedLat()    const { return standard_lon() && shifted_lat();  }

protected:

    bool global_lonlat() const {
      return domain().global()
        &&   not projection()
        &&   yspace().type() == "linear";
    }

    bool standard_lon() const {
        return x(0) == 0.;
    }

    bool standard_lat() const {
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

}
}
