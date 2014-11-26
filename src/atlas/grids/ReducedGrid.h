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

#include <eckit/memory/Builder.h>

#include "atlas/Grid.h"
#include "atlas/Util.h"

namespace atlas {
namespace grids {

/// @brief Reduced Grid
///
/// This class is a base class for all grids that can be described by
/// constant latitudes with a uniform distribution of points per latitude
/// in zonal direction.
/// This means any full grid and reduced grid, both regular, gaussian or other
/// such distribution can be represented with this class
class ReducedGrid: public Grid {
public:
  typedef eckit::SharedPtr<ReducedGrid> Ptr;

  static ReducedGrid* create( const eckit::Params& );
  static ReducedGrid* create( const GridSpec& );
  static ReducedGrid* create( const std::string& uid );

public:

  ReducedGrid();

  ReducedGrid( const eckit::Params& );

  ReducedGrid( const std::vector<size_t>& nlon, const std::vector<double>& lats );

  ReducedGrid( const int npts_per_lat[], const double lats[], const int nlat );

  static std::string className();

  virtual std::string uid() const;

  virtual std::string hash() const;

  virtual BoundBox bounding_box() const;

  virtual size_t npts() const;

  virtual void lonlat( double[] ) const;

  virtual void lonlat( std::vector<Point>& ) const;

  virtual std::string grid_type() const;

  virtual GridSpec spec() const;

  virtual bool same( const Grid& ) const;

  // number of latitudes in hemispherecd
  /// @todo temporary: this should go to the GaussianGrid classes
  int N() const { return N_; }

  int nlat() const;

  int nlon( int jlat ) const;

  int nlonmax() const;

  const std::vector<int>& npts_per_lat() const;

  const std::vector<double>& latitudes() const;

  double lon( const int jlat, const int jlon ) const;

  double lat( const int jlat ) const;

  void lonlat( const int jlon, const int jlat, double crd[] ) const;

  /// @brief Crop the grid according to the bounding box (or mask)
  void crop(const BoundBox&);
  void crop(const eckit::Params&);

protected:
  void setup(const eckit::Params&);
  void setup(const int nlat, const int npts_per_lat[], const double lats[]);
  void setup_colat_hemisphere(const int N, const int lon[], const double colat[], const AngleUnit);
  void setup_lat_hemisphere(const int N, const int lon[], const double lat[], const AngleUnit);

protected:
  std::vector<double> lat_;    ///<! Latitude values
  std::vector<int>    nlons_;  ///<! Number of points per latitude (int32 type for Fortran interoperability)
  std::vector<double> lonmin_; ///<! Value of minimum longitude per latitude [default=0]
  std::vector<double> lonmax_; ///<! Value of maximum longitude per latitude [default=~360 (one increment smaller)]
  size_t              npts_;   ///<! Total number of unique points in the grid
  int                 nlonmax_;
  BoundBox            bounding_box_;   ///<! bounding box for data, only points within are considered part of grid
  std::string         uid_;
  std::string         hash_;
  std::string         grid_type_;
  int                 N_;
};

register_BuilderT1(Grid,ReducedGrid,"reduced");


extern "C"
{
  int  atlas__ReducedGrid__nlat(ReducedGrid* This);
  void atlas__ReducedGrid__nlon(ReducedGrid* This, const int* &nlon, int &size);
  int atlas__ReducedGrid__npts(ReducedGrid* This);
  double atlas__ReducedGrid__lon(ReducedGrid* This,int jlat,int jlon);
  double atlas__ReducedGrid__lat(ReducedGrid* This,int jlat);
  void atlas__ReducedGrid__latitudes(ReducedGrid* This, const double* &lats, int &size);
}

} // namespace grids
} // namespace atlas

#endif // ReducedGrid_h
