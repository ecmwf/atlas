/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include <memory>

#include "eckit/memory/Builder.h"
#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/util/Config.h"

#include "atlas/grid/Spacing.h"


namespace atlas {
namespace grid {
namespace detail {
namespace grid {

/**
 * @brief Structured Grid
 *
 * This class is a base class for all grids that can be described by
 * constant latitudes with a uniform distribution of points per latitude
 * in zonal direction.
 * This means any full grid and reduced grid, both regular, gaussian or other
 * such distribution can be represented with this class
 */
class Structured : public Grid {

public:

  class Iterator: public Grid::Iterator {
  public:
    Iterator(const Structured& grid, bool begin = true):
        grid_(grid),
        i_(0),
        j_( begin ? 0 : grid.ny() ) {
    }

    virtual bool next(PointXY& xy) {

       if( j_<grid_.ny() && i_<grid_.nlon(j_) ) {

         xy = grid_.xy(j_,i_++);

         if( i_==grid_.nlon(j_) ) {
           j_++;
           i_=0;
         }
         return true;
       }
       return false;
    }
    
    
    virtual const PointXY operator *() const {
        return grid_.xy(j_,i_);
    }

    virtual const Grid::Iterator& operator ++() {
        ++i_;
        if( i_ == grid_.nlon(j_) ) {
          ++j_;
          i_=0;
        }
        return *this;
    }

    virtual bool operator ==(const Grid::Iterator &other) const {
        return j_ == static_cast<const Iterator&>(other).j_ && i_ == static_cast<const Iterator&>(other).i_;
    }

    virtual bool operator !=(const Grid::Iterator &other) const {
        return i_ != static_cast<const Iterator&>(other).i_ || j_ != static_cast<const Iterator&>(other).j_;
    }
    
    
  private:
    const Structured& grid_;
    size_t i_;
    size_t j_;
  };

public:

    struct XSpace {

      XSpace(long ny);

      size_t ny;

      // Minimum number of points across parallels (constant y)
      size_t nxmin;

      // Maximum number of points across parallels (constant y)
      size_t nxmax;

      /// Number of points per latitude
      std::vector<long> nx;

      /// Value of minimum longitude per latitude [default=0]
      std::vector<double> xmin;

      /// Value of maximum longitude per latitude [default=0]
      std::vector<double> xmax;

      /// Value of longitude increment
      std::vector<double> dx;
    };

    using YSpace = Spacing;

public:

    static std::string static_type();

public:

    Structured();

    Structured( Projection, XSpace*, YSpace, Domain );

    virtual ~Structured();

    virtual size_t npts() const {
        return npts_;
    }

    virtual void lonlat(std::vector<Point>&) const;

    //virtual std::string type() const {
    //    return grid_type_;
     //}

    virtual Spec spec() const;

    /**
     * Human readable name
     * @note: may not be unique, such as when reduced Gaussian grids have the same N numbers but different distribution of latitude points
     */
    virtual std::string name() const;

    virtual std::string type() const;

    virtual std::string getOptimalMeshGenerator() const {
        return "Structured";
    }

    virtual size_t N() const {
        return N_;
    }

    inline size_t ny() const {
        return y_.size();
    }

    inline size_t nlat() const {
        return y_.size();
    }

    inline size_t nlon( size_t jlat ) const {
        return static_cast<size_t>(nx_[jlat]);
    }

    inline size_t nlonmax() const {
        return nxmax_;
    }

    inline size_t nlonmin() const {
        return nxmin_;
    }

    inline const std::vector<long>& pl() const {
        return nx_;
    }

    inline const std::vector<double>& latitudes() const {
        return y_;
    }

    inline double lon( const size_t jlat, const size_t jlon ) const {
        return xmin_[jlat] + static_cast<double>(jlon) * dx_[jlat];
    }

    inline double lat( const size_t jlat ) const {
        return y_[jlat];
    }

    inline void lonlat( const size_t jlat, const size_t jlon, double crd[] ) const {
        crd[0] = lon(jlat,jlon);
        crd[1] = lat(jlat);
    }

    inline void xy( const size_t jlat, const size_t jlon, double p[] ) const {
        p[0] = lon(jlat,jlon);
        p[1] = lat(jlat);
    }

    PointXY xy( const size_t jlat, const size_t jlon ) const {
        return PointXY( lon(jlat,jlon), lat(jlat) );
    }

    PointLonLat geolonlat( const size_t jlon, const size_t jlat ) const {
        return projection_.lonlat( xy(jlat,jlon) );
    }

    void geoLonlat(const size_t jlon, const size_t jlat, PointLonLat &Pll) const {
      Pll.assign( projection_.lonlat( xy(jlat,jlon) ) );
    }

    inline bool reduced() const {
        return nlonmax() != nlonmin();
    }

    bool isPeriodicX() const { return periodic_x_; }

    const YSpace& yspace() const { return yspace_; }

    virtual Iterator* begin() const{ return new Iterator(*this); }
    virtual Iterator* end()   const{ return new Iterator(*this,false); }


protected: // methods

    virtual void print(std::ostream&) const;

    /// Hash of the PL array
    virtual void hash(eckit::MD5&) const;

    // void setup_cropped(const size_t ny, const double y[], const long nx[], const double xmin[], const double xmax[], const domain::Domain& dom);

    void compute_true_periodicity();

    void setup(
        const YSpace&              yspace,
        const std::vector<long>&   nx,
        const std::vector<double>& xmin,
        const std::vector<double>& xmax,
        const std::vector<double>& dx );

    void setup(
        const YSpace&     yspace,
        const long&       nx,
        const double&     xmin,
        const double&     xmax,
        const double&     dx );

protected:

    /// Number of latitudes in hemisphere (only makes sense for global grids)
    size_t N_;

    // Minimum number of points across parallels (constant y)
    size_t nxmin_;

    // Maximum number of points across parallels (constant y)
    size_t nxmax_;

    /// Total number of unique points in the grid
    size_t npts_;

    /// Latitude values
    std::vector<double> y_;

    /// Number of points per latitude
    std::vector<long> nx_;

    /// Value of minimum longitude per latitude [default=0]
    std::vector<double> xmin_;

    /// Value of maximum longitude per latitude [default=0]
    std::vector<double> xmax_;

    /// Value of longitude increment
    std::vector<double> dx_;

    /// Periodicity in x-direction
    bool periodic_x_;

private:
    std::unique_ptr<XSpace> xspace_;
    YSpace yspace_;
    mutable std::string type_;
};


extern "C"
{
    void atlas__grid__Structured__delete(Structured* This);
    const Structured *atlas__grid__Structured(char* identifier);
    const Structured* atlas__grid__Structured__config(util::Config* conf);
    Structured* atlas__grid__CustomStructured_int(size_t nlat, double lat[], int nlon[]);
    Structured* atlas__grid__CustomStructured_long(size_t nlat, double lat[], long nlon[]);
    Structured* atlas__grid__CustomStructured_lonmin_lonmax_int(size_t nlat, double lat[], int nlon[], double lonmin[], double lonmax[]);
    Structured* atlas__grid__CustomStructured_lonmin_lonmax_long(size_t nlat, double lat[], long nlon[], double lonmin[], double lonmax[]);
    Structured* atlas__grid__regular__RegularGaussian(size_t N);
    Structured* atlas__grid__reduced__ReducedGaussian_int(size_t N, int nlon[]);
    Structured* atlas__grid__reduced__ReducedGaussian_long(size_t N, long nlon[]);
    Structured* atlas__grid__regular__RegularLonLat(size_t nlon, size_t nlat);
    Structured* atlas__grid__regular__ShiftedLonLat(size_t nlon, size_t nlat);
    Structured* atlas__grid__regular__ShiftedLon(size_t nlon, size_t nlat);
    Structured* atlas__grid__regular__ShiftedLat(size_t nlon, size_t nlat);

    void   atlas__grid__Structured__pl        (Structured* This, const long* &pl, size_t &size);
    size_t atlas__grid__Structured__N         (Structured* This);
    size_t atlas__grid__Structured__nlat      (Structured* This);
    size_t atlas__grid__Structured__nlon      (Structured* This, size_t jlat);
    size_t atlas__grid__Structured__nlonmin   (Structured* This);
    size_t atlas__grid__Structured__nlonmax   (Structured* This);
    size_t atlas__grid__Structured__npts      (Structured* This);
    double atlas__grid__Structured__lat       (Structured* This, size_t jlat);
    double atlas__grid__Structured__lon       (Structured* This, size_t jlat, size_t jlon);
    void   atlas__grid__Structured__lonlat    (Structured* This, size_t jlat, size_t jlon, double crd[]);
    void   atlas__grid__Structured__latitudes (Structured* This, const double* &lats, size_t &size);
    int    atlas__grid__Structured__reduced   (Structured* This);
}


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
