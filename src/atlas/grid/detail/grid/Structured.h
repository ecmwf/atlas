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

  class IteratorXY: public Grid::IteratorXY {
  public:
    IteratorXY(const Structured& grid, bool begin = true):
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

    virtual const Grid::IteratorXY& operator ++() {
        ++i_;
        if( i_ == grid_.nx(j_) ) {
          ++j_;
          i_=0;
        }
        return *this;
    }

    virtual bool operator ==(const Grid::IteratorXY &other) const {
        return j_ == static_cast<const IteratorXY&>(other).j_ && i_ == static_cast<const IteratorXY&>(other).i_;
    }

    virtual bool operator !=(const Grid::IteratorXY &other) const {
        return i_ != static_cast<const IteratorXY&>(other).i_ || j_ != static_cast<const IteratorXY&>(other).j_;
    }


  private:
    const Structured& grid_;
    size_t i_;
    size_t j_;
  };

  class IteratorLonLat: public Grid::IteratorLonLat {
  public:
    IteratorLonLat(const Structured& grid, bool begin = true):
        grid_(grid),
        i_(0),
        j_( begin ? 0 : grid.ny() ) {
    }

    virtual bool next(PointLonLat& lonlat) {

       if( j_<grid_.ny() && i_<grid_.nx(j_) ) {

         lonlat = grid_.lonlat(i_++,j_);

         if( i_==grid_.nx(j_) ) {
           j_++;
           i_=0;
         }
         return true;
       }
       return false;
    }


    virtual const PointLonLat operator *() const {
        return grid_.lonlat(i_,j_);
    }

    virtual const Grid::IteratorLonLat& operator ++() {
        ++i_;
        if( i_ == grid_.nx(j_) ) {
          ++j_;
          i_=0;
        }
        return *this;
    }

    virtual bool operator ==(const Grid::IteratorLonLat &other) const {
        return j_ == static_cast<const IteratorLonLat&>(other).j_ && i_ == static_cast<const IteratorLonLat&>(other).i_;
    }

    virtual bool operator !=(const Grid::IteratorLonLat &other) const {
        return i_ != static_cast<const IteratorLonLat&>(other).i_ || j_ != static_cast<const IteratorLonLat&>(other).j_;
    }


  private:
    const Structured& grid_;
    size_t i_;
    size_t j_;
  };

public:

    class XSpace {

        class Implementation : public eckit::Owned {

        public:

            Implementation( const std::array<double,2>& interval, const std::vector<long>& N, bool endpoint=true );

            Implementation( const Spacing& );

            Implementation( const Config& );

            Implementation( const std::vector<Config>& );

            size_t ny() const { return ny_; }

            // Minimum number of points across parallels (constant y)
            size_t nxmin() const { return nxmin_; }

            // Maximum number of points across parallels (constant y)
            size_t nxmax() const { return nxmax_; }

            /// Number of points per latitude
            const std::vector<long>& nx() const { return nx_; }

            /// Value of minimum longitude per latitude [default=0]
            const std::vector<double>& xmin() const { return xmin_; }

            /// Value of maximum longitude per latitude [default=0]
            const std::vector<double>& xmax() const { return xmax_; }

            /// Value of longitude increment
            const std::vector<double>& dx() const { return dx_; }

            Spec spec() const;

        private:

            void reserve( long ny );

        private:

            size_t ny_;
            size_t nxmin_;
            size_t nxmax_;
            std::vector<long> nx_;
            std::vector<double> xmin_;
            std::vector<double> xmax_;
            std::vector<double> dx_;
        };

    public:

        XSpace();

        XSpace( const XSpace& );

        XSpace( const Spacing& );

        XSpace( const std::array<double,2>& interval, const std::vector<long>& N, bool endpoint=true );

        XSpace( const Config& );

        XSpace( const std::vector<Config>& );

        size_t ny() const { return impl_->ny(); }

        // Minimum number of points across parallels (constant y)
        size_t nxmin() const { return impl_->nxmin(); }

        // Maximum number of points across parallels (constant y)
        size_t nxmax() const { return impl_->nxmax(); }

        /// Number of points per latitude
        const std::vector<long>& nx() const { return impl_->nx(); }

        /// Value of minimum longitude per latitude [default=0]
        const std::vector<double>& xmin() const { return impl_->xmin(); }

        /// Value of maximum longitude per latitude [default=0]
        const std::vector<double>& xmax() const { return impl_->xmax(); }

        /// Value of longitude increment
        const std::vector<double>& dx() const { return impl_->dx(); }

        Spec spec() const { return impl_->spec(); }

    private:

        eckit::SharedPtr<Implementation> impl_;
    };

    using YSpace = Spacing;

public:

    static std::string static_type();

public:

    Structured( const std::string&, XSpace, YSpace, Projection, Domain );
    Structured( XSpace, YSpace, Projection, Domain );

    virtual ~Structured();

    virtual size_t size() const {
        return npts_;
    }

    virtual Spec spec() const;

    /**
     * Human readable name
     * @note: may not be unique, such as when reduced Gaussian grids have the same N numbers but different distribution of latitude points
     */
    virtual std::string name() const;

    virtual std::string type() const;

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

    inline double x( size_t i, size_t j ) const {
        return xmin_[j] + static_cast<double>(i) * dx_[j];
    }

    inline double y( size_t j ) const {
        return y_[j];
    }

    inline void xy( size_t i, size_t j, double crd[] ) const {
        crd[0] = x(i,j);
        crd[1] = y(j);
    }

    PointXY xy( size_t i, size_t j ) const {
        return PointXY( x(i,j), y(j) );
    }

    PointLonLat lonlat( size_t i, size_t j ) const {
        return projection_.lonlat( xy(i,j) );
    }

    void lonlat(size_t i, size_t j, double crd[]) const {
      xy(i,j,crd);
      projection_.xy2lonlat(crd);
    }

    inline bool reduced() const {
        return nxmax() != nxmin();
    }

    bool periodic() const { return periodic_x_; }

    const XSpace& xspace() const { return xspace_; }
    const YSpace& yspace() const { return yspace_; }

    virtual IteratorXY* xy_begin() const{ return new IteratorXY(*this); }
    virtual IteratorXY* xy_end()   const{ return new IteratorXY(*this,false); }
    virtual IteratorLonLat* lonlat_begin() const{ return new IteratorLonLat(*this); }
    virtual IteratorLonLat* lonlat_end()   const{ return new IteratorLonLat(*this,false); }


protected: // methods

    virtual void print(std::ostream&) const;

    /// Hash of the PL array
    virtual void hash(eckit::MD5&) const;

    void computeTruePeriodicity();

    void computeDomain();

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
    XSpace xspace_;
    YSpace yspace_;
    mutable std::string type_;
};


extern "C"
{
    void atlas__grid__Structured__delete(Structured* This);
    const Structured *atlas__grid__Structured(char* identifier);
    const Structured* atlas__grid__Structured__config(util::Config* conf);
    Structured* atlas__grid__regular__RegularGaussian(long N);
    Structured* atlas__grid__reduced__ReducedGaussian_int(long N, int nlon[]);
    Structured* atlas__grid__reduced__ReducedGaussian_long(long N, long nlon[]);
    Structured* atlas__grid__regular__RegularLonLat(long nlon, long nlat);
    Structured* atlas__grid__regular__ShiftedLonLat(long nlon, long nlat);
    Structured* atlas__grid__regular__ShiftedLon(long nlon, long nlat);
    Structured* atlas__grid__regular__ShiftedLat(long nlon, long nlat);

    void   atlas__grid__Structured__nx_array  (Structured* This, const long* &pl, size_t &size);
    long   atlas__grid__Structured__nx        (Structured* This, long j);
    long   atlas__grid__Structured__ny        (Structured* This);
    long   atlas__grid__Structured__nxmin     (Structured* This);
    long   atlas__grid__Structured__nxmax     (Structured* This);
    long   atlas__grid__Structured__size      (Structured* This);
    double atlas__grid__Structured__y         (Structured* This, long j);
    double atlas__grid__Structured__x         (Structured* This, long i, long j);
    void   atlas__grid__Structured__xy        (Structured* This, long i, long j, double crd[]);
    void   atlas__grid__Structured__lonlat    (Structured* This, long i, long j, double crd[]);
    void   atlas__grid__Structured__y_array   (Structured* This, const double* &lats, size_t &size);
    int    atlas__grid__Structured__reduced   (Structured* This);

    long   atlas__grid__Gaussian__N           (Structured* This);

}


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
