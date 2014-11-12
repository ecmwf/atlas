/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef ReducedGrid_h
#define ReducedGrid_h

#include "atlas/Grid.h"

namespace atlas {

/// @brief Reduced Grid
///
/// This class is a base class for all grids that can be described by
/// constant latitudes with a uniform distribution of points on each latitude
/// in zonal direction.
/// This means any full grid and reduced grid, both regular, gaussian or other
/// distribution can be represented with this class
class ReducedGrid: public Grid {
public:

  ReducedGrid();

  ReducedGrid( const std::vector<size_t>& nlon, const std::vector<double>& lats );

  ReducedGrid( const int nlons[], const double lats[], const int nlat );

  static std::string className();

  virtual std::string uid() const;

  virtual std::string hash() const;

  virtual BoundBox boundingBox() const;

  virtual size_t nPoints() const;

  virtual void coordinates( std::vector<double>& ) const;

  virtual void coordinates( std::vector<Point>& ) const;

  virtual std::string gridType() const;

  virtual GridSpec spec() const;

  virtual bool same( const Grid& ) const;

  int nlat() const;

  int nlon( int jlat ) const;

  int nlonmax() const;

  const std::vector<int>& nlons() const;

  const std::vector<double>& lat() const;

  double lon( const int jlat, const int jlon ) const;

  double lat( const int jlat ) const;

  void lonlat( const int jlon, const int jlat, double crd[] ) const;

protected:
  void setup(const int nlat, const int nlons[], const double lats[]);

protected:
  std::vector<double> lat_;   ///<! Latitude values
  std::vector<int>    nlons_;  ///<! Number of points per latitude (int32 type for Fortran interoperability)
  size_t              npts_;   ///<! Total number of unique points in the grid
  int                 nlonmax_;
  BoundBox            bbox_;   ///<! bounding box for data, only points within are considered part of grid
};

} // namespace atlas

#endif // ReducedGrid_h
