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
    ComputePointXY( const CubedSphere& grid ) : grid_( grid ) {}
    void operator()( idx_t i, idx_t j, idx_t t, PointXY& point ) {
      grid_.xy( i, j, t, point.data() );
    }
    const CubedSphere& grid_;
  };

  struct ComputePointLonLat {
    ComputePointLonLat( const CubedSphere& grid ) : grid_( grid ) {}
    void operator()( idx_t i, idx_t j, idx_t t, PointLonLat& point ) {
      grid_.lonlat( i, j, t, point.data() );
    }
    const CubedSphere& grid_;
  };

  // -----------------------------------------------------------------------------------------------

  template <typename Base, typename ComputePoint>
  class CubedSphereIterator : public Base {
  public:
    CubedSphereIterator( const CubedSphere& grid, bool begin = true ) :
                          grid_( grid ), i_( begin ? 0 : grid_.CubeNx() ),
                          j_( begin ? 0 : grid_.CubeNx() ), t_( begin ? 0 : 5 ),
                          compute_point{grid_} {
      if ( grid_.inGrid(i_, j_, t_) ) {
        compute_point( i_, j_, t_, point_ );
      }
    }

    virtual bool next( typename Base::value_type& point ) {
      if ( grid_.inGrid(i_, j_, t_) && !grid_.finalElement(i_, j_, t_)) {
        compute_point( i_, j_, t_, point );
        int* ijt = grid_.nextElement(i_, j_, t_);
        i_ = ijt[0];
        j_ = ijt[1];
        t_ = ijt[2];
        return true;
      }
      return false;
    }

    virtual const typename Base::reference operator*() const { return point_; }

    virtual const Base& operator++() {
      int* ijt = grid_.nextElement(i_, j_, t_);
      i_ = ijt[0];
      j_ = ijt[1];
      t_ = ijt[2];
      compute_point( i_, j_, t_, point_ );
      return *this;
    }

    virtual const Base& operator+=( typename Base::difference_type distance ) {
      idx_t d = distance;
      for (int n = 0; n < d; n++) {
        int* ijt = grid_.nextElement(i_, j_, t_);
        i_ = ijt[0];
        j_ = ijt[1];
        t_ = ijt[2];
      }
      compute_point( i_, j_, t_, point_ );
      return *this;
    }

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
        int* ijt = grid_.nextElement(i, j, t);
        i = ijt[0];
        j = ijt[1];
        t = ijt[2];
        ++d;
      }
      ATLAS_ASSERT( !found, "CubedSphereIterator.distance: cycled entire grid without finding other" );
      return d;
    }

    virtual bool operator==( const Base& other ) const {
      return i_ == static_cast<const CubedSphereIterator&>( other ).i_ &&
             j_ == static_cast<const CubedSphereIterator&>( other ).j_ &&
             t_ == static_cast<const CubedSphereIterator&>( other ).t_;
    }

    virtual bool operator!=( const Base& other ) const {
      return i_ != static_cast<const CubedSphereIterator&>( other ).i_ ||
             j_ != static_cast<const CubedSphereIterator&>( other ).j_ ||
             t_ != static_cast<const CubedSphereIterator&>( other ).t_;
    }

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
  using IteratorXY     = CubedSphereIterator<Grid::IteratorXY, ComputePointXY>;
  using IteratorLonLat = CubedSphereIterator<Grid::IteratorLonLat, ComputePointLonLat>;

  static std::string static_type();

  // Constructors
  CubedSphere( const std::string&, int, Projection );
  CubedSphere( int, Projection );
  CubedSphere( const CubedSphere& );

  // Functions overloaded from the base class
  virtual ~CubedSphere() override;
  virtual idx_t size() const override {
    return accumulate(npts_.begin(), npts_.end(), 0);
  }
  virtual Spec spec() const override;
  virtual std::string name() const override;
  virtual std::string type() const override;

  // Return number of faces on cube
  inline idx_t GetCubeNx() const { return CubeNx_; }
  inline idx_t CubeNx() const { return CubeNx_; }


  // Tole specific access to x and y locations
  // -----------------------------------------

  inline double x123( idx_t i, idx_t j, idx_t t ) const {
    return static_cast<double>(xs_[t]) + static_cast<double>(i);
  }

  inline double x456( idx_t i, idx_t j, idx_t t ) const {
    return static_cast<double>(xs_[t])  + (y123( i, j, t )-static_cast<double>(ys_[t]));
  }

  inline double y123( idx_t i, idx_t j, idx_t t ) const {
    return static_cast<double>(ys_[t]) + static_cast<double>(j);
  }

  inline double y456( idx_t i, idx_t j, idx_t t ) const {
    return static_cast<double>(ysr_[t]) - x123( i, j, t ) + static_cast<double>(xs_[t]);
  }

  // Lambdas for access to appropriate functions for tile
  // ----------------------------------------------------

  std::vector<std::function<double(int, int, int)>> xtile =
    {[this](int i, int j, int t) -> double{return this->x123(i, j, t);},
     [this](int i, int j, int t) -> double{return this->x456(i, j, t);}
    };

  std::vector<std::function<double(int, int, int)>> ytile =
    {[this](int i, int j, int t) -> double{return this->y123(i, j, t);},
     [this](int i, int j, int t) -> double{return this->y456(i, j, t);}
    };

  // Functions for returning xy
  // --------------------------

  inline void xy( idx_t i, idx_t j, idx_t t, double xyt[] ) const {
    xyt[0] = xtile.at(round(t/6.0))(i, j, t);
    xyt[1] = ytile.at(round(t/6.0))(i, j, t);
    xyt[2] = static_cast<double>(t);
  }

  PointXY xy( idx_t i, idx_t j, idx_t t ) const {
    return PointXY( xtile.at(round(t/6.0))(i, j, t), ytile.at(round(t/6.0))(i, j, t) );
  }

  // Functions for returning lonlat
  // ------------------------------

  void lonlat( idx_t i, idx_t j, idx_t t, double lonlat[] ) const {
    double xytll[5];
    xytll[0] = i;
    xytll[1] = j;
    xytll[2] = t;
    projection_.xy2lonlat( xytll );
    lonlat[LON] = xytll[3+LON];
    lonlat[LAT] = xytll[3+LAT];
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

  int* nextElement(const idx_t i, const idx_t j, const idx_t t) const {

    int* ijt = new int[100];

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

    // Internal points
    if (ijt[0] == CubeNx_-1 && ijt[1] < CubeNx_-1 && ijt[2] == 6) {  // Final point
      ijt[0] = 0;
      ijt[1] = 0;
      ijt[2] = 0;
    } else if (ijt[0] == CubeNx_-1 && ijt[1] == CubeNx_-1) {  // Corner
      ijt[0] = 0;
      ijt[1] = 0;
      ijt[2] = ijt[2] + 1;
    } else if (ijt[0] == CubeNx_-1) {  // Edge
      ijt[0] = 0;
      ijt[1] = ijt[1] + 1;
    } else {
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

  // Nuber of faces on tile
  idx_t CubeNx_;

  // Start points in x,y direction
  int xs_[6];
  int ys_[6];
  int ysr_[6]; // (for panels 4, 5, 6)

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
