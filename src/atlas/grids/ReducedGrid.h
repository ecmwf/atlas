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
#include "eckit/config/Parametrisation.h"
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

  /// FIXME: ReducedGrid should not be instantiatable.
  ///        Only leaf classes should be instantiatable.
  ///        This constructor should be used only by derived types
  ReducedGrid(const Domain& d = Domain::makeGlobal());

  ReducedGrid( const eckit::Parametrisation& );

  ReducedGrid( const std::vector<double>& lats,
               const std::vector<size_t>& nlon,
               const Domain& d = Domain::makeGlobal());

  ReducedGrid( size_t nlat,
               const double lats[],
               const int npts_per_lat[],
               const Domain& d = Domain::makeGlobal());

  virtual BoundBox boundingBox() const;

  virtual size_t npts() const;

  virtual void lonlat( std::vector<Point>& ) const;

  virtual std::string gridType() const;

  virtual GridSpec spec() const;

  /// number of latitudes in hemisphere
  virtual size_t N() const;

  size_t nlat() const;

  size_t nlon( size_t jlat ) const;

  size_t nlonmax() const;

  const std::vector<int>& npts_per_lat() const;

  const std::vector<double>& latitudes() const;

  double lon( const size_t jlat, const size_t jlon ) const;

  double lat( const size_t jlat ) const;

  void lonlat( const size_t jlon, const size_t jlat, double crd[] ) const;

private: // methods

  virtual size_t copyLonLatMemory(double* pts, size_t size) const;

  virtual void print(std::ostream&) const;

  /// Human readable name
  /// May not be unique, especially when reduced gauss grids have the same N numbers
  /// but different distribution of latitude points
  virtual std::string shortName() const;

protected:

  /// Hash of the PL array
  virtual void hash(eckit::MD5&) const;

  /// @note Domain is already set when calling setup()
  void setup(const eckit::Parametrisation& );
  /// @note Domain is already set when calling setup()
  void setup( const size_t nlat, const double lats[], const int npts_per_lat[] );
  /// @note Domain is already set when calling setup()
  void setup( const size_t nlat, const double lats[], const int nlons[], const double lonmin[], const double lonmax[] );
  /// @note Domain is already set when calling setup()
  void setup_lat_hemisphere( const size_t N, const double lat[], const int lon[], const AngleUnit );

protected:

  size_t              N_;
  size_t              nlonmax_;

  size_t              npts_;          ///<! Total number of unique points in the grid

  std::string         grid_type_;

  std::string         shortName_;

  std::vector<double> lat_;    ///<! Latitude values
  std::vector<int>    nlons_;  ///<! Number of points per latitude (int32 type for Fortran interoperability)
  std::vector<double> lonmin_; ///<! Value of minimum longitude per latitude [default=0]
  std::vector<double> lonmax_; ///<! Value of maximum longitude per latitude [default=~360 (one increment smaller)]

  BoundBox            bounding_box_;  ///<! bounding box cache

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
