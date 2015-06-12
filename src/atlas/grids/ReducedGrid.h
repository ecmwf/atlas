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

#include "eckit/memory/Builder.h"

#include "atlas/Parameters.h"
#include "atlas/Grid.h"

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

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
  static ReducedGrid* create( const std::string& shortName );

public:

  static std::string className();
  static std::string grid_type_str() { return "reduced"; }

  ReducedGrid();

  ReducedGrid( const eckit::Params& );

  ReducedGrid( const std::vector<double>& lats, const std::vector<size_t>& nlon );

  ReducedGrid( const int nlat, const double lats[], const int npts_per_lat[] );

  virtual BoundBox bounding_box() const;

  virtual size_t npts() const;

  virtual void lonlat( double[] ) const;

  virtual void lonlat( std::vector<Point>& ) const;

  virtual std::string grid_type() const;

  virtual GridSpec spec() const;

  /// number of latitudes in hemisphere
  virtual int N() const;

  int nlat() const;

  int nlon( int jlat ) const;

  int nlonmax() const;

  const std::vector<int>& npts_per_lat() const;

  const std::vector<double>& latitudes() const;

  double lon( const int jlat, const int jlon ) const;

  double lat( const int jlat ) const;

  void lonlat( const int jlon, const int jlat, double crd[] ) const;

  /// @brief Mask the grid according to the domain
  virtual void mask( const Domain& );
  virtual void mask( const eckit::Params& );

private: // methods

  /// Human readable name
  /// May not be unique, especially when reduced gauss grids have the same N numbers
  /// but different distribution of latitude points
  virtual std::string shortName() const;

protected:

  /// Hash of the PL array
  virtual void hash(eckit::MD5&) const;

  void setup( const eckit::Params& );
  void setup( const int nlat, const double lats[], const int npts_per_lat[] );
  void setup( const int nlat, const double lats[], const int nlons[], const double lonmin[], const double lonmax[] );
  void setup_lat_hemisphere( const int N, const double lat[], const int lon[], const AngleUnit );

protected:

  int                 N_;
  int                 nlonmax_;

  size_t              npts_;          ///<! Total number of unique points in the grid

  BoundBox            bounding_box_;  ///<! bounding box for data, only points within are considered part of grid

  std::string         grid_type_;

  std::string         shortName_;

  std::vector<double> lat_;    ///<! Latitude values
  std::vector<int>    nlons_;  ///<! Number of points per latitude (int32 type for Fortran interoperability)
  std::vector<double> lonmin_; ///<! Value of minimum longitude per latitude [default=0]
  std::vector<double> lonmax_; ///<! Value of maximum longitude per latitude [default=~360 (one increment smaller)]

};

//------------------------------------------------------------------------------------------------------

extern "C"
{
  int  atlas__ReducedGrid__nlat(ReducedGrid* This);
  void atlas__ReducedGrid__nlon(ReducedGrid* This, const int* &nlon, int &size);
  int atlas__ReducedGrid__npts(ReducedGrid* This);
  double atlas__ReducedGrid__lon(ReducedGrid* This,int jlat,int jlon);
  double atlas__ReducedGrid__lat(ReducedGrid* This,int jlat);
  void atlas__ReducedGrid__latitudes(ReducedGrid* This, const double* &lats, int &size);
}

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

#endif // ReducedGrid_h
