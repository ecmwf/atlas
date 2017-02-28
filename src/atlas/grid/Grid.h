#pragma once

#include <string>

#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/grid/detail/grid/Regular.h"
#include "atlas/grid/Projection.h"
#include "atlas/grid/Domain.h"

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------

class Grid {

public:

  using Config = atlas::util::Config;
  using grid_t = detail::grid::Grid;
  using Spec = eckit::Properties;

public:

    Grid();
    Grid( const Grid& );
    Grid( const grid_t* );
    Grid( const std::string& );
    Grid( const Config& );

    operator bool() const { return grid_; }
    operator const detail::grid::Grid&() const { return *grid_.get(); }
    bool operator==( const Grid& other ) const { return grid_ == other.grid_; }
    bool operator!=( const Grid& other ) const { return grid_ != other.grid_; }

    size_t npts() const { return grid_->npts(); }

    const Projection& projection() const { return grid_->projection(); }
    const Domain& domain() const { return grid_->domain(); }
    std::string name() { return grid_->shortName(); }
    std::string uid() { return grid_->uniqueId(); }

    Spec spec() const { return grid_->spec(); }

    // TODO: replace with iterator
    void lonlat(std::vector<PointLonLat>& points) const { return grid_->lonlat(points); }
    void fillLonLat(std::vector<double>& vector) const { return grid_->fillLonLat(vector); }
    void fillLonLat(double array[], size_t arraySize) const { return grid_->fillLonLat(array,arraySize); }

private:

    eckit::SharedPtr<const grid_t> grid_;

protected:

    const grid_t* raw() { return grid_.get(); }
};

//---------------------------------------------------------------------------------------------------------------------

class StructuredGrid: public Grid {

public:

  using structured_t = detail::grid::Structured;
  using YSpace = structured_t::YSpace;

public:

    StructuredGrid();
    StructuredGrid( const Grid& );
    StructuredGrid( const grid_t* );
    StructuredGrid( const std::string& );
    StructuredGrid( const Config& );

    operator bool() const { return grid_; }

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

    bool periodic() const { return grid_->isPeriodicX(); }

    const YSpace& yspace() const { return grid_->yspace(); }

private:

    static structured_t* create( const Config& );
    const structured_t* grid_;
};

//---------------------------------------------------------------------------------------------------------------------

class RegularGrid: public StructuredGrid {

public:

    using regular_t = structured_t;

public:

    RegularGrid( const Grid& );
    RegularGrid( const grid_t* );
    RegularGrid( const std::string& );
    RegularGrid( const Config& );

    operator bool() const { return grid_; }

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

    static regular_t* create( const Config& );
    const regular_t* grid_;
    size_t nx_ = {0};
};

//---------------------------------------------------------------------------------------------------------------------

}
}
