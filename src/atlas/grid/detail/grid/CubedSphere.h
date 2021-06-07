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
#include "atlas/runtime/Log.h"
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
      grid_.lonlat( i, j, t, point.data() );
    }
    const CubedSphere& grid_;
  };

  class PointTIJ : public std::array<idx_t,3> {
  public:
      using std::array<idx_t,3>::array;
      idx_t t() const { return data()[0]; }
      idx_t i() const { return data()[1]; }
      idx_t j() const { return data()[2]; }
      idx_t& t() { return data()[0]; }
      idx_t& i() { return data()[1]; }
      idx_t& j() { return data()[2]; }
  };

  struct ComputePointTIJ {
      ComputePointTIJ( const CubedSphere& grid ) : grid_( grid ) {}
      void operator()( idx_t i, idx_t j, idx_t t, PointTIJ& point ) {
          point.t() = t;
          point.i() = i;
          point.j() = j;
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
                          grid_( grid ), i_( begin ? 0 : grid_.N() ),
                          j_( begin ? 0 : grid_.N() ), t_( begin ? 0 : 5 ),
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

  class IteratorTIJ_Base : public IteratorT<IteratorTIJ_Base, PointTIJ> {};
  using IteratorTIJ = CubedSphereIterator<IteratorTIJ_Base, ComputePointTIJ>;

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
  inline idx_t GetN() const { return N_; }
  inline idx_t N() const { return N_; }

  // Return number of tiles
  inline idx_t GetNTiles() const { return nTiles_; }

  void xy2xyt(const double xy[], double xyt[]) const;

  void xyt2xy(const double xyt[], double xy[]) const;

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

  std::vector<std::function<double(int, int, int)>> xtile;
  std::vector<std::function<double(int, int, int)>> ytile;

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
    this->xy(i, j, t, lonlat);  // outputing xy in lonlat
    projection_.xy2lonlat( lonlat ); // converting xy to lonlat
  }

  PointLonLat lonlat( idx_t i, idx_t j, idx_t t ) const {
    double lonlat[2];
    this->lonlat(i, j, t, lonlat);
    return PointLonLat( lonlat[LON], lonlat[LAT] );
  }

  // Check on whether i, j, t values are for extra point on tile 1
  // -------------------------------------------------------------

  inline bool extraPoint1(idx_t i, idx_t j, idx_t t) const {
    if (i == 0 && j == N_ && t == 0) {
      return true;
    } else {
      return false;
    }
  }

  // Check on whether i, j, t values are for extra point on tile 2
  // -------------------------------------------------------------

  inline bool extraPoint2(idx_t i, idx_t j, idx_t t) const {
    if (i == N_ && j == 0 && t == 1) {
      return true;
    } else {
      return false;
    }
  }

  // Check whether i, j, t is in grid
  // --------------------------------

  inline bool inGrid(idx_t i, idx_t j, idx_t t) const {
    if (t >= 0 && t <= 5) {
      if (i >= 0 && i <= N_-1) {
        if (i >= 0 && i <= N_-1) {
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
    if (i == N_-1 && j == N_-1 && t == 5) {
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
    if (ijt[2] == 0 && ijt[0] == N_-1 && ijt[1] == N_-1) {
      ijt[0] = 0;
      ijt[1] = N_;
      return ijt;
    }

    // Moving from extra point 2
    if (extraPoint2(ijt[0], ijt[1], ijt[2])) {
      ijt[0] = 0;
      ijt[1] = 1;
      return ijt;
    }

    // Moving to extra point 2
    if (ijt[0] == N_-1 && ijt[1] == 0 && ijt[2] == 1) {
      ijt[0] = N_;
      return ijt;
    }


    if (ijt[0] == N_-1 && ijt[1] == N_-1 && ijt[2] == nTiles_-1) {  // Final point
      ijt[0] = N_;
      ijt[1] = N_;
      ijt[2] = nTiles_-1;
    } else if (ijt[0] == N_-1 && ijt[1] == N_-1) {  // Corner
      ijt[0] = 0;
      ijt[1] = 0;
      ijt[2] = ijt[2] + 1;
    } else if (ijt[0] == N_-1) {  // Edge
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
  virtual std::unique_ptr<IteratorTIJ> tij_begin() const {
      return std::unique_ptr<IteratorTIJ>( new IteratorTIJ( *this ) );
  }
  virtual std::unique_ptr<IteratorTIJ> tij_end() const {
      return std::unique_ptr<IteratorTIJ>( new IteratorTIJ( *this, false ) );
  }

  // Default configurations

  Config meshgenerator() const override;

  Config partitioner() const override;

protected:
  virtual void print( std::ostream& ) const override;

  virtual void hash( eckit::Hash& ) const override;

  virtual RectangularLonLatDomain lonlatBoundingBox() const override;

  Domain computeDomain() const;

  // Number of faces on tile
  idx_t N_;

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
};


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
