/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>
#include <memory>
#include <iostream>

#include "atlas/array.h"
#include "atlas/grid/Spacing.h"
#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

/**
 * @brief CubedSphere Grid
 *
 * This class is a base class for all grids that can be described as
 * a cubed sphere.
 * This means the equidistant, equiangular or any other
 * such distribution can be represented with this class
 */
class CubedSphere : public Grid {
private:
    struct ComputePointXY {
        ComputePointXY( const CubedSphere& grid ) : grid_( grid ), ny_( grid_.ny() ) {}
        void operator()( idx_t i, idx_t j, idx_t t, PointXY& point ) {
            if ( j < ny_ ) {  // likely
                grid_.xy( i, j, t, point.data() );
            }
        }
        const CubedSphere& grid_;
        idx_t ny_;
    };
    struct ComputePointLonLat {
        ComputePointLonLat( const CubedSphere& grid ) : grid_( grid ), ny_( grid_.ny() ) {}
        void operator()( idx_t i, idx_t j, idx_t t, PointLonLat& point ) {
            if ( j < ny_ ) {  // likely
                grid_.lonlat( i, j, t, point.data() );
            }
        }
        const CubedSphere& grid_;
        idx_t ny_;
    };

    template <typename Base, typename ComputePoint>
    class CubedSphereIterator : public Base {
    public:
        CubedSphereIterator( const CubedSphere& grid, bool begin = true ) :
            grid_( grid ),
            ny_( grid_.ny() ),
            i_( 0 ),
            j_( begin ? 0 : grid_.ny() ),
            t_( 0 ),
            compute_point{grid_} {
            if ( j_ != ny_ && grid_.size() ) {
                compute_point( i_, j_, t_, point_ );
            }
        }
        virtual bool next( typename Base::value_type& point ) {
            if ( j_ < ny_ && i_ < grid_.nx( j_ ) ) {
                compute_point( i_++, j_, t_, point );

                if ( i_ == grid_.nx( j_ ) ) {
                    j_++;
                    i_ = 0;
                }
                return true;
            }
            return false;
        }

        virtual const typename Base::reference operator*() const { return point_; }

        virtual const Base& operator++() {
            ++i_;
            if ( i_ == grid_.nx( j_ ) ) {
                ++j_;
                i_ = 0;
            }
            compute_point( i_, j_, t_, point_ );
            return *this;
        }

        virtual const Base& operator+=( typename Base::difference_type distance ) {
            idx_t d = distance;
            while ( j_ != ny_ && d >= ( grid_.nx( j_ ) - i_ ) ) {
                d -= ( grid_.nx( j_ ) - i_ );
                ++j_;
                i_ = 0;
            }
            i_ += d;
            compute_point( i_, j_, t_, point_ );
            return *this;
        }

        virtual typename Base::difference_type distance( const Base& other ) const {
            const auto& _other               = static_cast<const CubedSphereIterator&>( other );
            typename Base::difference_type d = 0;
            idx_t j                          = j_;
            idx_t i                          = i_;
            while ( j < _other.j_ ) {
                d += grid_.nx( j ) - i;
                ++j;
                i = 0;
            }
            d += _other.i_;
            return d;
        }

        virtual bool operator==( const Base& other ) const {
            return j_ == static_cast<const CubedSphereIterator&>( other ).j_ &&
                   i_ == static_cast<const CubedSphereIterator&>( other ).i_;
        }

        virtual bool operator!=( const Base& other ) const {
            return i_ != static_cast<const CubedSphereIterator&>( other ).i_ ||
                   j_ != static_cast<const CubedSphereIterator&>( other ).j_;
        }

        virtual std::unique_ptr<Base> clone() const {
            auto result    = new CubedSphereIterator( grid_, false );
            result->i_     = i_;
            result->j_     = j_;
            result->point_ = point_;
            return std::unique_ptr<Base>( result );
        }

    public:
        const CubedSphere& grid_;
        idx_t ny_;
        idx_t i_;
        idx_t j_;
        idx_t t_;
        typename Base::value_type point_;
        ComputePoint compute_point;
    };

public:
    using IteratorXY     = CubedSphereIterator<Grid::IteratorXY, ComputePointXY>;
    using IteratorLonLat = CubedSphereIterator<Grid::IteratorLonLat, ComputePointLonLat>;

public:

public:
    static std::string static_type();

public:
    CubedSphere( const std::string&, int, Projection );
    CubedSphere( int, Projection );
    CubedSphere( const CubedSphere& );

    virtual ~CubedSphere() override;

    virtual idx_t size() const override { return npts_; }

    virtual Spec spec() const override;

    /**
   * Human readable name
   * Either the name is the one given at construction as a canonical named grid,
   * or the name "cubedsphere"
   */
    virtual std::string name() const override;

    virtual std::string type() const override;

    inline idx_t ny() const { return static_cast<idx_t>( y_.size() ); }

    inline idx_t nx( idx_t j ) const { return static_cast<idx_t>( nx_[j] ); }

    inline const std::vector<idx_t>& nx() const { return nx_; }

    inline const std::vector<double>& y() const { return y_; }

    inline double dx( idx_t j ) const { return 1.0; }

    inline double x( idx_t i, idx_t j, idx_t t ) const {
      int pos = xs_[t] + i;
      double posd = static_cast<double>(pos);
      return posd;
    }

    inline double y( idx_t i, idx_t j, idx_t t ) const {
      int pos = ys_[t] + j;
      double posd = static_cast<double>(pos);
      return posd;
    }
    inline void xy( idx_t i, idx_t j, idx_t t, double crd[] ) const {
      crd[0] = x( i, j, t );
      crd[1] = y( i, j, t );
      crd[2] = static_cast<double>(t);
    }

    PointXY xy( idx_t i, idx_t j, idx_t t ) const { return PointXY( x( i, j, t ), y( i, j, t ) ); }

    PointLonLat lonlat( idx_t i, idx_t j, idx_t t ) const { return projection_.lonlat( xy( i, j, t ) ); }

    void lonlat( idx_t i, idx_t j, idx_t t, double crd[] ) const {
      double xyt[3];
      xy( i, j, t, xyt );

      double xytll[5];
      xytll[0] = i;
      xytll[1] = j;
      xytll[2] = t;

      projection_.xy2lonlat( xytll );
    }

    virtual std::unique_ptr<Grid::IteratorXY> xy_begin() const override {
        return std::unique_ptr<Grid::IteratorXY>( new IteratorXY( *this ) );
    }
    virtual std::unique_ptr<Grid::IteratorXY> xy_end() const override {
        return std::unique_ptr<Grid::IteratorXY>( new IteratorXY( *this, false ) );
    }
    virtual std::unique_ptr<Grid::IteratorLonLat> lonlat_begin() const override {
        return std::unique_ptr<Grid::IteratorLonLat>( new IteratorLonLat( *this ) );
    }
    virtual std::unique_ptr<Grid::IteratorLonLat> lonlat_end() const override {
        return std::unique_ptr<Grid::IteratorLonLat>( new IteratorLonLat( *this, false ) );
    }

    int GetCubeNx() const { return CubeNx_; }

protected:  // methods
    virtual void print( std::ostream& ) const override;

    virtual void hash( eckit::Hash& ) const override;

    virtual RectangularLonLatDomain lonlatBoundingBox() const override;

    Domain computeDomain() const;

protected:
    // Total number of unique points in the grid
    idx_t npts_;

    // Latitude values
    std::vector<double> y_;

    // Number of points per latitude
    std::vector<idx_t> nx_;

    // Value of minimum longitude per latitude [default=0]
    std::vector<double> xmin_;

    // Value of maximum longitude per latitude [default=0]
    std::vector<double> xmax_;

    // Nuber of faces on tile
    int CubeNx_;

    // Start points in x,y direction
    int xs_[6];
    int ys_[6];

    // End points in x,y direction
    int xe_[6];
    int ye_[6];

    array::ArrayT<int> xArray_();

private:
    std::string name_ = {"cubedsphere"};
    mutable std::string type_;
};

extern "C" {
void atlas__grid__CubedSphere__delete( CubedSphere* This );
const CubedSphere* atlas__grid__CubedSphere( char* identifier );
const CubedSphere* atlas__grid__CubedSphere__config( util::Config* conf );


void atlas__grid__CubedSphere__nx_array( CubedSphere* This, const idx_t*& nx, idx_t& size );
idx_t atlas__grid__CubedSphere__nx( CubedSphere* This, idx_t j );
idx_t atlas__grid__CubedSphere__ny( CubedSphere* This );
idx_t atlas__grid__CubedSphere__size( CubedSphere* This );
double atlas__grid__CubedSphere__y( CubedSphere* This, idx_t i, idx_t j, idx_t t );
double atlas__grid__CubedSphere__x( CubedSphere* This, idx_t i, idx_t j, idx_t t );
void atlas__grid__CubedSphere__xy( CubedSphere* This, idx_t i, idx_t j, idx_t t, double crd[] );
void atlas__grid__CubedSphere__lonlat( CubedSphere* This, idx_t i, idx_t j, idx_t t, double crd[] );
void atlas__grid__CubedSphere__y_array( CubedSphere* This, const double*& lats, idx_t& size );

//idx_t atlas__grid__Gaussian__N( CubedSphere* This );
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
