#pragma once

#include <string>

#include "atlas/grid/Grid.h"
#include "atlas/grid/Structured.h"
#include "atlas/grid/Regular.h"
#include "atlas/grid/ptr/Projection.h"
#include "atlas/grid/ptr/Domain.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {
namespace ptr {

//---------------------------------------------------------------------------------------------------------------------

class Grid {

public:

    Grid();
    Grid( const Grid& );
    Grid( const atlas::grid::Grid* );
    Grid( const std::string& );

    operator bool() const { return grid_; }
    operator const atlas::grid::Grid&() const { return *grid_.get(); }

    const Projection& projection() { return projection_; }
    const Domain& domain() { return domain_; }

private:

    eckit::SharedPtr<const atlas::grid::Grid> grid_;
    Projection projection_;
    Domain domain_;

protected:

    const atlas::grid::Grid* raw() { return grid_.get(); }
};

//---------------------------------------------------------------------------------------------------------------------

class Structured: public Grid {

public:

    Structured();
    Structured( const Grid& );
    Structured( const atlas::grid::Grid* );
    Structured( const eckit::Parametrisation& );

    operator bool() const { return grid_; }
    //operator const atlas::grid::Grid*() const { return grid_; }
    operator Grid() const { return Grid(grid_); }

    inline size_t ny() const {
        return grid_->ny();
    }

    inline size_t nx( size_t j ) const {
        return grid_->nlon(j);;
    }

    inline const std::vector<long>& nx() const {
        return grid_->pl();
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

    PointXY xy( const size_t i, const size_t j ) const {
      return PointXY( x(i,j), y(j ) );
    }

    PointLonLat lonlat( const size_t i, const size_t j ) const {
      return grid_->geolonlat(i,j);
    }

    inline bool reduced() const {
        return grid_->reduced();
    }

    bool periodic() const { return grid_->isPeriodicX(); }

private:

    static atlas::grid::Structured* create( const eckit::Parametrisation& );
    const atlas::grid::Structured* grid_;
};

//---------------------------------------------------------------------------------------------------------------------

class Regular: public Structured {

public:

    Regular( const Grid& );
    Regular( const atlas::grid::Grid* );
    Regular( const eckit::Parametrisation& );

    operator bool() const { return grid_; }

    size_t nx() const { return nx_; }

    inline double x( const size_t i ) const {
        return Structured::x(i,0);
    }

    using Structured::y;

    inline double y( const size_t j ) const {
        return grid_->lat(j);
    }

    PointXY xy( const size_t i, const size_t j ) const {
      return PointXY( x(i), y(j ) );
    }

private:

    static atlas::grid::Structured* create( const eckit::Parametrisation& );
    const atlas::grid::Structured* grid_;
    size_t nx_ = {0};
};

//---------------------------------------------------------------------------------------------------------------------

}
}
}
