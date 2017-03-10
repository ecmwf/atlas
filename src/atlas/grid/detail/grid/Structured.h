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
#include <array>

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

       if( j_<grid_.ny() && i_<grid_.nx(j_) ) {

         xy = grid_.xy(i_++,j_);

         if( i_==grid_.nx(j_) ) {
           j_++;
           i_=0;
         }
         return true;
       }
       return false;
    }


    virtual const PointXY operator *() const {
        return grid_.xy(i_,j_);
    }

    virtual const Grid::Iterator& operator ++() {
        ++i_;
        if( i_ == grid_.nx(j_) ) {
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

      XSpace( const std::array<double,2>& interval, const std::vector<long>& N, bool endpoint=true );

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

    Structured( const std::string&, Projection, XSpace*, YSpace, Domain );
    Structured( Projection, XSpace*, YSpace, Domain );

    virtual ~Structured();

    virtual size_t npts() const {
        return npts_;
    }

    virtual Spec spec() const;

    /**
     * Human readable name
     * @note: may not be unique, such as when reduced Gaussian grids have the same N numbers but different distribution of latitude points
     */
    virtual std::string name() const;

    virtual std::string type() const;

    virtual std::string getOptimalMeshGenerator() const {
        return "structured";
    }

    inline size_t ny() const {
        return y_.size();
    }

    inline size_t nx( size_t j ) const {
        return static_cast<size_t>(nx_[j]);
    }

    inline size_t nxmax() const {
        return nxmax_;
    }

    inline size_t nxmin() const {
        return nxmin_;
    }

    inline const std::vector<long>& nx() const {
        return nx_;
    }

    inline const std::vector<double>& y() const {
        return y_;
    }

    inline double x( const size_t i, const size_t j ) const {
        return xmin_[j] + static_cast<double>(i) * dx_[j];
    }

    inline double y( const size_t j ) const {
        return y_[j];
    }

    inline void xy( const size_t i, const size_t j, double crd[] ) const {
        crd[0] = x(i,j);
        crd[1] = y(j);
    }

    PointXY xy( const size_t i, const size_t j ) const {
        return PointXY( x(i,j), y(j) );
    }

    PointLonLat lonlat( const size_t i, const size_t j ) const {
        return projection_.lonlat( xy(i,j) );
    }

    void lonlat(const size_t i, const size_t j, PointLonLat &Pll) const {
      Pll.assign( projection_.lonlat( xy(i,j) ) );
    }

    inline bool reduced() const {
        return nxmax() != nxmin();
    }

    bool periodic() const { return periodic_x_; }

    const YSpace& yspace() const { return yspace_; }

    virtual Iterator* begin() const{ return new Iterator(*this); }
    virtual Iterator* end()   const{ return new Iterator(*this,false); }


protected: // methods

    virtual void print(std::ostream&) const;

    /// Hash of the PL array
    virtual void hash(eckit::MD5&) const;

    void compute_true_periodicity();

    void crop( const Domain& );

protected:

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
    std::string name_ = {"structured"};
    std::unique_ptr<XSpace> xspace_;
    YSpace yspace_;
    mutable std::string type_;
};


extern "C"
{
    void atlas__grid__Structured__delete(Structured* This);
    const Structured *atlas__grid__Structured(char* identifier);
    const Structured* atlas__grid__Structured__config(util::Config* conf);
    Structured* atlas__grid__CustomStructured_int(long nlat, double lat[], int nlon[]);
    Structured* atlas__grid__CustomStructured_long(long nlat, double lat[], long nlon[]);
    Structured* atlas__grid__CustomStructured_lonmin_lonmax_int(long nlat, double lat[], int nlon[], double lonmin[], double lonmax[]);
    Structured* atlas__grid__CustomStructured_lonmin_lonmax_long(long nlat, double lat[], long nlon[], double lonmin[], double lonmax[]);
    Structured* atlas__grid__regular__RegularGaussian(long N);
    Structured* atlas__grid__reduced__ReducedGaussian_int(long N, int nlon[]);
    Structured* atlas__grid__reduced__ReducedGaussian_long(long N, long nlon[]);
    Structured* atlas__grid__regular__RegularLonLat(long nlon, long nlat);
    Structured* atlas__grid__regular__ShiftedLonLat(long nlon, long nlat);
    Structured* atlas__grid__regular__ShiftedLon(long nlon, long nlat);
    Structured* atlas__grid__regular__ShiftedLat(long nlon, long nlat);

    void   atlas__grid__Structured__pl        (Structured* This, const long* &pl, size_t &size);
    long   atlas__grid__Structured__nlat      (Structured* This);
    long   atlas__grid__Structured__nlon      (Structured* This, long jlat);
    long   atlas__grid__Structured__nlonmin   (Structured* This);
    long   atlas__grid__Structured__nlonmax   (Structured* This);
    long   atlas__grid__Structured__npts      (Structured* This);
    double atlas__grid__Structured__lat       (Structured* This, long jlat);
    double atlas__grid__Structured__lon       (Structured* This, long jlat, long jlon);
    void   atlas__grid__Structured__lonlat    (Structured* This, long jlat, long jlon, double crd[]);
    void   atlas__grid__Structured__latitudes (Structured* This, const double* &lats, size_t &size);
    int    atlas__grid__Structured__reduced   (Structured* This);

    long   atlas__grid__Gaussian__N           (Structured* This);

}


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
