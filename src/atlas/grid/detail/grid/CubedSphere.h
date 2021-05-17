/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#include "atlas/array.h"
#include "atlas/grid/Spacing.h"
#include "atlas/grid/detail/grid/Grid.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"
#include "atlas/util/Point.h"

#include "atlas/projection/detail/CubedSphereProjectionBase.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

/**
 * @brief CubedSphere Grid
 *
 * This class is a base class for all grids that can be described as
 * a cubed sphere.
 *
 * For more detail on this implementation see atlas/grid/CubedSphereGrid.h
 */

class CubedSphere : public Grid {
private:

  // Get the position in the xy plane and return as PointXY object
  struct ComputePointXY {
    ComputePointXY( const CubedSphere& grid ) : grid_( grid ) {}
    void operator()( idx_t i, idx_t j, idx_t t, PointXY& point ) {
      grid_.xy( i, j, t, point.data() );
    }
    const CubedSphere& grid_;
  };

  // Get the lonlat and return as PointLonLat object
  struct ComputePointLonLat {
    ComputePointLonLat( const CubedSphere& grid ) : grid_( grid ) {}
    void operator()( idx_t i, idx_t j, idx_t t, PointLonLat& point ) {
      std::cout << "i,j,t,point.data" << i << " " << j << " " << t << " " << point << std::endl;
      grid_.lonlat( i, j, t, point.data() );
    }
    const CubedSphere& grid_;
  };

  // -----------------------------------------------------------------------------------------------

  template <typename Base, typename ComputePoint>
  class CubedSphereIterator : public Base {
  public:

    // Create an iterator and return the first or last point. If begin is true it starts at the
    // beginning of the iterator, otherwise at the end. Class is templated and point can be xy or
    // lonlat.
    CubedSphereIterator( const CubedSphere& grid, bool begin = true ) :
                          grid_( grid ), i_( begin ? 0 : grid_.CubeNx() ),
                          j_( begin ? 0 : grid_.CubeNx() ), t_( begin ? 0 : 5 ),
                          compute_point{grid_} {
      // Check that point lies in grid and if so return the xy/lonlat
      if ( grid_.inGrid(i_, j_, t_) ) {
        compute_point( i_, j_, t_, point_ );
      }
    }

    // Return the point and move iterator to the next location
    virtual bool next( typename Base::value_type& point ) {
      if ( grid_.inGrid(i_, j_, t_) && !grid_.finalElement(i_, j_, t_)) {
        compute_point( i_, j_, t_, point );
        std::unique_ptr<int[]> ijt = grid_.nextElement(i_, j_, t_);
        i_ = ijt[0];
        j_ = ijt[1];
        t_ = ijt[2];
        return true;
      }
      return false;
    }

    // * operator
    virtual const typename Base::reference operator*() const { return point_; }

    // ++ operator, move to next element in grid iterator and return point
    virtual const Base& operator++() {
      std::unique_ptr<int[]> ijt = grid_.nextElement(i_, j_, t_);
      i_ = ijt[0];
      j_ = ijt[1];
      t_ = ijt[2];
      compute_point( i_, j_, t_, point_ );
      return *this;
    }

    // += operator, move some distance d through the iterator and return point
    virtual const Base& operator+=( typename Base::difference_type distance ) {
      idx_t d = distance;
      for (int n = 0; n < d; n++) {
       std::unique_ptr<int[]> ijt = grid_.nextElement(i_, j_, t_);
        i_ = ijt[0];
        j_ = ijt[1];
        t_ = ijt[2];
      }
      compute_point( i_, j_, t_, point_ );
      return *this;
    }

    // Given two positions in the grid iterator return the distance, which for the cubed-sphere
    // grid is just the number of grid points between the two points.
    virtual typename Base::difference_type distance( const Base& other ) const {
      const auto& _other = static_cast<const CubedSphereIterator&>( other );
      typename Base::difference_type d = 0;
      idx_t i = i_;
      idx_t j = j_;
      idx_t t = t_;
      bool found = false;
      for (int n = 0; n < grid_.size(); n++) {
        if (i == _other.i_ && j == _other.j_ && t == _other.t_) {
          found = true;
          break;
        }
        std::unique_ptr<int[]> ijt = grid_.nextElement(i, j, t);
        i = ijt[0];
        j = ijt[1];
        t = ijt[2];
        ++d;
      }
      ATLAS_ASSERT( !found, "CubedSphereIterator.distance: cycled entire grid without finding other" );
      return d;
    }

    // == operator for checking two positions in the iterator are equal
    virtual bool operator==( const Base& other ) const {
      return i_ == static_cast<const CubedSphereIterator&>( other ).i_ &&
             j_ == static_cast<const CubedSphereIterator&>( other ).j_ &&
             t_ == static_cast<const CubedSphereIterator&>( other ).t_;
    }

    // != operator for checking that two positions in the iterator are not equal
    virtual bool operator!=( const Base& other ) const {
      return i_ != static_cast<const CubedSphereIterator&>( other ).i_ ||
             j_ != static_cast<const CubedSphereIterator&>( other ).j_ ||
             t_ != static_cast<const CubedSphereIterator&>( other ).t_;
    }

    // Clone the grid iterator
    virtual std::unique_ptr<Base> clone() const {
      auto result    = new CubedSphereIterator( grid_, false );
      result->i_     = i_;
      result->j_     = j_;
      result->t_     = t_;
      result->point_ = point_;
      return std::unique_ptr<Base>( result );
    }

    const CubedSphere& grid_;
    idx_t i_;
    idx_t j_;
    idx_t t_;
    typename Base::value_type point_;
    ComputePoint compute_point;
  };

  // -----------------------------------------------------------------------------------------------

public:

  // Iterators for returning xy or lonlat
  using IteratorXY     = CubedSphereIterator<Grid::IteratorXY, ComputePointXY>;
  using IteratorLonLat = CubedSphereIterator<Grid::IteratorLonLat, ComputePointLonLat>;

  static std::string static_type();

  // Constructors
  CubedSphere( const std::string&, int, Projection );
  CubedSphere( int, Projection );
  CubedSphere( const CubedSphere& );

  // Destructor
  virtual ~CubedSphere() override;

  // Return total grid size
  virtual idx_t size() const override {
    return accumulate(npts_.begin(), npts_.end(), 0);
  }

  // Return information about the grid
  virtual Spec spec() const override;
  virtual std::string name() const override;
  virtual std::string type() const override;

  // Return number of faces on cube
  inline idx_t GetCubeNx() const { return CubeNx_; }
  inline idx_t CubeNx() const { return CubeNx_; }

  // Return number of tiles
  inline idx_t GetNTiles() const { return nTiles_; }

  inline idx_t tile(const double xy[] ) const {
    // Assume one face-edge is of length 90 degrees.
    //
    //   y ^
    //     |
    //    135           -------
    //     |           |       |
    //     |           |   2   |
    //     |           |       |
    //     45   ------- ------- ------- -------
    //     |   |       |       |       |       |
    //     |   |   0   |   1   |   3   |   4   |
    //     |   |       |       |       |       |
    //    -45   -------  ------ ------- -------
    //     |                           |       |
    //     |                           |   5   |
    //     |                           |       |
    //   -135                           -------
    //     ----0-------90------180-----270----360--->  x

    /*
    const double x = xy[0]/90.;
    const double& y = xy[1];
    if( x < 2. ) {
      return  y > 45. ? 2 : std::floor(x);
    }
    else {
      return y < -45. ? 5 : std::floor(x+1.);
    }
    */

    idx_t t{-1};

    if ((xy[0] >= 0.) && ( xy[1] >= -45.) && (xy[0] < 90.) && (xy[1] < 45.)) {
       t = 0;
    } else if ((xy[0] >= 90.) && ( xy[1] >= -45.) && (xy[0] < 180.) && (xy[1] < 45.)) {
       t = 1;
    } else if ((xy[0] >= 90.) && ( xy[1] >= 45.) && (xy[0] < 180.) && (xy[1] < 135.)) {
       t = 2;
    } else if ((xy[0] >= 180.) && ( xy[1] > -45.) && (xy[0] < 270.) && (xy[1] <= 45.)) {
       t = 3;
    } else if ((xy[0] >= 270.) && ( xy[1] > -45.) && (xy[0] < 360.) && (xy[1] <= 45.)) {
       t = 4;
    } else if ((xy[0] >= 270.) && ( xy[1] > -135.) && (xy[0] < 360.) && (xy[1] <= -45.)) {
       t = 5;
    }

    // extra points
    if ((xy[0] == 0.) && (xy[1] == 45.)) t = 0;
    if ((xy[0] == 180.) && (xy[1] == -45.)) t = 1;

    return t;

  }

  void xy2xyt(const double xy[], double xyt[]) const {
      // xy is in degrees while xyt is in radians
      // (alpha, beta) and tiles.

      double normalisedX = xy[0]/90.;
      double normalisedY = (xy[1] + 135.)/90.;

      double CubeNxDouble = static_cast<double>(CubeNx_);

      std::vector<double> yOffset{CubeNxDouble,
                                  CubeNxDouble,
                                  2. *  CubeNxDouble,
                                  CubeNxDouble,
                                  CubeNxDouble,
                                  0};


      //  xyt[0] = (normalisedX - std::floor(normalisedX))* M_PI_2 - M_PI_4;
      //  xyt[1] = (normalisedY - std::floor(normalisedY))* M_PI_2 - M_PI_4;
      xyt[0] =  (normalisedX - std::floor(normalisedX)) * static_cast<double>(CubeNx_)
            +  xs_[static_cast<size_t>(xyt[2])];

      xyt[1] = (normalisedY - std::floor(normalisedY)) * static_cast<double>(CubeNx_)
            + yOffset[static_cast<size_t>(xyt[2])];

      xyt[2] = tile(xy);
  }

  void xyt2xy(const double xyt[], double xy[]) const {
      // xy is in degrees
      // while xyt is in number of grid points
      // (alpha, beta) and tiles.
      std::vector<double> xOffsetDeg{0., 90., 90., 180, 270, 270};
      std::vector<double> yOffsetDeg{-45., -45, 45, -45, -45, -135};


      double N = static_cast<double>(CubeNx_);
      std::vector<double> xOffsetIndex{0, N, N, 2*N, 3*N,  3*N};
      std::vector<double> yOffsetIndex{N, N, 2*N, N,  N, 0};

      double normalisedX =
       (xyt[0] -  xOffsetIndex[static_cast<size_t>(xyt[2])])/N;
      double normalisedY =
       (xyt[1] - yOffsetIndex[static_cast<size_t>(xyt[2])])/N;
      xy[0] = normalisedX * 90. + xOffsetDeg[xyt[2]];
      xy[1] = normalisedY * 90. + yOffsetDeg[xyt[2]];
      std::cout << "xyt xy xs_ ysr_ :: "  << xyt[0] << " " << xyt[1] << " " << xyt[2]
                << " " << xy[0] << " " << xy[1] << " " << xs_[static_cast<size_t>(xyt[2])]
                << " " << ysr_[static_cast<size_t>(xyt[2])]  << " "
                << ys_[static_cast<size_t>(xyt[2])]   << std::endl;
  }

  // Tile specific access to x and y locations
  // -----------------------------------------

  inline double x123( idx_t i, idx_t t ) const {
    return static_cast<double>(xs_[t]) + static_cast<double>(i);
  }

  inline double x456( idx_t j, idx_t t ) const {
    return static_cast<double>(xs_[t])  + (y123( j, t )-static_cast<double>(ys_[t]));
  }

  inline double y123( idx_t j, idx_t t ) const {
    return static_cast<double>(ys_[t]) + static_cast<double>(j);
  }

  inline double y456( idx_t i, idx_t t ) const {
    return static_cast<double>(ysr_[t]) - x123( i, t ) + static_cast<double>(xs_[t]);
  }

  // Lambdas for access to appropriate functions for tile
  // ----------------------------------------------------

  std::vector<std::function<double(int, int, int)>> xtile =
    {[this](int i, int j, int t) {return this->x123(i, t);},
     [this](int i, int j, int t) {return this->x456(j, t);}
    };

  std::vector<std::function<double(int, int, int)>> ytile =
    {[this](int i, int j, int t) {return this->y123(j, t);},
     [this](int i, int j, int t) {return this->y456(i, t);}
    };

  // Functions for returning xy
  // --------------------------

  inline void xyt( idx_t i, idx_t j, idx_t t, double crd[] ) const {
    std::size_t tIndex = static_cast<std::size_t>(tileCases_ * t / nTiles_);
    crd[0] = xtile.at(tIndex)(i, j, t);
    crd[1] = ytile.at(tIndex)(i, j, t);
    crd[2] = static_cast<double>(t);
  }

  PointXY xyt( idx_t i, idx_t j, idx_t t ) const {
    std::size_t tIndex = static_cast<std::size_t>(tileCases_ * t / nTiles_);
    return PointXY( xtile.at(tIndex)(i, j, t), ytile.at(tIndex)(i, j, t) );
  }

  inline void xy( idx_t i, idx_t j, idx_t t, double xy[] ) const {
    double crd[3];
    this->xyt(i,j,t,crd);
    std::cout << "   " << std::endl;
    std::cout <<"crd " << crd[0] << " " << crd[1]  << " " << crd[2] << std::endl;
    this->xyt2xy(crd, xy);
  }

  PointXY xy( idx_t i, idx_t j, idx_t t ) const {
    double crd[2];
    this->xy(i,j,t,crd);
    return PointXY(crd[0], crd[1]);
  }

  // Functions for returning lonlat, either as array or PointLonLat
  // --------------------------------------------------------------

  void lonlat( idx_t i, idx_t j, idx_t t, double lonlat[] ) const {

    this->xy(i, j, t, lonlat);  // outputing xy in lonlat C array
    std::cout << "detail::CubedSphere before xy2lonlat :: ijt xy = " << i << " " << j << " " << t << " " << lonlat[0] << " " << lonlat[1] << std::endl;
    projection_.xy2lonlat( lonlat ); // converting xy to lonlat
    std::cout << "detail::CubedSphere after xy2lonlat :: ijt lonlat = " << i << " " << j << " " << t << " " <<  lonlat[0] << " " << lonlat[1] << std::endl;

  }

  PointLonLat lonlat( idx_t i, idx_t j, idx_t t ) const {
    double lonlat[2];
    this->lonlat(i, j, t, lonlat);
    return PointLonLat( lonlat[LON], lonlat[LAT] );
  }


  // Check on whether i, j, t values are for extra point on tile 1
  // -------------------------------------------------------------

  inline bool extraPoint1(idx_t i, idx_t j, idx_t t) const {
    if (i == 0 && j == CubeNx_ && t == 0) {
      return true;
    } else {
      return false;
    }
  }

  // Check on whether i, j, t values are for extra point on tile 2
  // -------------------------------------------------------------

  inline bool extraPoint2(idx_t i, idx_t j, idx_t t) const {
    if (i == CubeNx_ && j == 0 && t == 1) {
      return true;
    } else {
      return false;
    }
  }

  // Check whether i, j, t is in grid
  // --------------------------------

  inline bool inGrid(idx_t i, idx_t j, idx_t t) const {
    if (t >= 0 && t <= 5) {
      if (i >= 0 && i <= CubeNx_-1) {
        if (i >= 0 && i <= CubeNx_-1) {
          return true;
        }
      }
    }
    if (extraPoint1(i, j, t)) {return true;}
    if (extraPoint2(i, j, t)) {return true;}
    return false;
  }

  // Check on whether the final element
  // ----------------------------------

  bool finalElement(idx_t i, idx_t j, idx_t t) const {
    if (i == CubeNx_-1 && j == CubeNx_-1 && t == 5) {
      return true;
    }
    return false;
  }

  // Move to next grid element in an iterator
  // ----------------------------------------

  std::unique_ptr<int[]> nextElement(const idx_t i, const idx_t j, const idx_t t) const {

    std::unique_ptr<int[]> ijt (new int[3]);

    ijt[0] = i;
    ijt[1] = j;
    ijt[2] = t;

    // Moving from extra point 1
    if (extraPoint1(ijt[0], ijt[1], ijt[2])) {
      ijt[0] = 0;
      ijt[1] = 0;
      ijt[2] = 1;
      return ijt;
    }

    // Moving to extra point 1
    if (ijt[2] == 0 && ijt[0] == CubeNx_-1 && ijt[1] == CubeNx_-1) {
      ijt[0] = 0;
      ijt[1] = CubeNx_;
      return ijt;
    }

    // Moving from extra point 2
    if (extraPoint2(ijt[0], ijt[1], ijt[2])) {
      ijt[0] = 0;
      ijt[1] = 1;
      return ijt;
    }

    // Moving to extra point 2
    if (ijt[0] == CubeNx_-1 && ijt[1] == 0 && ijt[2] == 1) {
      ijt[0] = CubeNx_;
      return ijt;
    }


    if (ijt[0] == CubeNx_-1 && ijt[1] == CubeNx_-1 && ijt[2] == nTiles_-1) {  // Final point
      ijt[0] = CubeNx_;
      ijt[1] = CubeNx_;
      ijt[2] = nTiles_-1;
    } else if (ijt[0] == CubeNx_-1 && ijt[1] == CubeNx_-1) {  // Corner
      ijt[0] = 0;
      ijt[1] = 0;
      ijt[2] = ijt[2] + 1;
    } else if (ijt[0] == CubeNx_-1) {  // Edge
      ijt[0] = 0;
      ijt[1] = ijt[1] + 1;
    } else { // Internal points
      ijt[0] = ijt[0] + 1;
    }

    return ijt;
  }

  // Iterator start/end positions
  // ----------------------------

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

protected:
  virtual void print( std::ostream& ) const override;

  virtual void hash( eckit::Hash& ) const override;

  virtual RectangularLonLatDomain lonlatBoundingBox() const override;

  Domain computeDomain() const;

  // Number of faces on tile
  idx_t CubeNx_;

  // Number of tiles
  static const idx_t nTiles_ = 6;

  const idx_t tileCases_ = 2;

  // Start points in x,y direction
  int xs_[nTiles_];
  int ys_[nTiles_];
  int ysr_[nTiles_]; // (for panels 4, 5, 6)

  // Number of unique points on each tile
  std::vector<int> npts_;

private:
  std::string name_ = {"cubedsphere"};
  mutable std::string type_;
};

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
