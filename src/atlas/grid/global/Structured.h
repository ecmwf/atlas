/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef Structured_h
#define Structured_h

#include "eckit/memory/Builder.h"
#include "eckit/config/Parametrisation.h"
#include "atlas/internals/Parameters.h"
#include "atlas/grid/global/Global.h"

namespace atlas {
namespace grid {
namespace global {

//------------------------------------------------------------------------------

/// @brief Reduced Grid
///
/// This class is a base class for all grids that can be described by
/// constant latitudes with a uniform distribution of points per latitude
/// in zonal direction.
/// This means any full grid and reduced grid, both regular, gaussian or other
/// such distribution can be represented with this class

class Structured: public Global {
public:

  typedef eckit::SharedPtr<Structured> Ptr;

  static Structured* create( const eckit::Parametrisation& );
  static Structured* create( const eckit::Properties& );
  static Structured* create( const std::string& shortName );

public:

  static std::string className();
  static std::string grid_type_str() { return "structured"; }

  Structured( const Domain& dom=Domain::makeGlobal() );

  virtual size_t npts() const;

  virtual void lonlat( std::vector<Point>& ) const;

  virtual std::string gridType() const;

  virtual eckit::Properties spec() const = 0;

  /// number of latitudes in hemisphere
  virtual size_t N() const;

  size_t nlat() const;

  size_t nlon( size_t jlat ) const;

  size_t nlonmax() const;

  size_t nlonmin() const;

  const std::vector<long>& pl() const;

  const std::vector<double>& latitudes() const;

  double lon( const size_t jlat, const size_t jlon ) const;

  double lat( const size_t jlat ) const;

  void lonlat( const size_t jlat, const size_t jlon, double crd[] ) const;

  bool reduced() const { return nlonmax() != nlonmin(); }

  /// Human readable name
  /// May not be unique, especially when reduced gauss grids have the same N numbers
  /// but different distribution of latitude points
  virtual std::string shortName() const;

  virtual std::string getOptimalMeshGenerator() const;

protected: // methods

  virtual size_t copyLonLatMemory(double* pts, size_t size) const;

  virtual void print(std::ostream&) const;

  /// Hash of the PL array
  virtual void hash(eckit::MD5&) const;

  void setup( const size_t nlat, const double lats[], const long pl[] );
  void setup( const size_t nlat, const double lats[], const long pl[], const double lonmin[] );
  void setup_lat_hemisphere( const size_t N, const double lat[], const long lon[] );

protected:

  size_t              N_;
  size_t              nlonmin_;
  size_t              nlonmax_;

  size_t              npts_;          ///<! Total number of unique points in the grid

  std::string         grid_type_; ///< to be removed -- only instantiate leaf classees

  std::string         shortName_; ///< to be removed -- only instantiate leaf classees

  std::vector<double> lat_;    ///<! Latitude values
  std::vector<long>   pl_;     ///<! Number of points per latitude

  std::vector<double> lonmin_; ///<! Value of minimum longitude per latitude [default=0]

private:
  std::vector<double> lon_inc_; ///<! Value of longitude increment
};

//------------------------------------------------------------------------------

inline size_t Structured::nlat() const
{
  return lat_.size();
}

inline size_t Structured::nlon(size_t jlat) const
{
  return pl_[jlat];
}

inline size_t Structured::nlonmin() const
{
  return nlonmin_;
}

inline size_t Structured::nlonmax() const
{
  return nlonmax_;
}

inline const std::vector<long>&  Structured::pl() const
{
  return pl_;
}

inline double Structured::lon(const size_t jlat, const size_t jlon) const
{
  return lonmin_[jlat] + (double)jlon * lon_inc_[jlat];;
}

inline double Structured::lat(const size_t jlat) const
{
  return lat_[jlat];
}

inline void Structured::lonlat( const size_t jlat, const size_t jlon, double crd[] ) const
{
  crd[0] = lon(jlat,jlon);
  crd[1] = lat(jlat);
}


//------------------------------------------------------------------------------
#define CONFIG eckit::Parametrisation
extern "C"
{
  void atlas__grid__global__Structured__delete(Structured* This);
  Structured* atlas__grid__global__Structured(char* identifier);
  Structured* atlas__grid__global__Structured__config(CONFIG* conf);
  Structured* atlas__grid__global__CustomStructured_int(size_t nlat, double lat[], int nlon[]);
  Structured* atlas__grid__global__CustomStructured_long(size_t nlat, double lat[], long nlon[]);
  Structured* atlas__grid__global__CustomStructured_lonmin_int(size_t nlat, double lat[], int nlon[], double lonmin[]);
  Structured* atlas__grid__global__CustomStructured_lonmin_long(size_t nlat, double lat[], long nlon[], double lonmin[]);
  Structured* atlas__grid__global__gaussian__RegularGaussian(size_t N);
  Structured* atlas__grid__global__gaussian__ReducedGaussian_int(size_t N, int nlon[]);
  Structured* atlas__grid__global__gaussian__ReducedGaussian_long(size_t N, long nlon[]);
  Structured* atlas__grid__global__lonlat__RegularLonLat(size_t nlon, size_t nlat);
  Structured* atlas__grid__global__lonlat__ShiftedLonLat(size_t nlon, size_t nlat);
  Structured* atlas__grid__global__lonlat__ShiftedLon(size_t nlon, size_t nlat);
  Structured* atlas__grid__global__lonlat__ShiftedLat(size_t nlon, size_t nlat);

  void      atlas__grid__global__Structured__pl       (Structured* This, const long* &pl, size_t &size);
  size_t    atlas__grid__global__Structured__N        (Structured* This);
  size_t    atlas__grid__global__Structured__nlat     (Structured* This);
  size_t    atlas__grid__global__Structured__nlon     (Structured* This, size_t jlat);
  size_t    atlas__grid__global__Structured__nlonmin  (Structured* This);
  size_t    atlas__grid__global__Structured__nlonmax  (Structured* This);
  size_t    atlas__grid__global__Structured__npts     (Structured* This);
  double atlas__grid__global__Structured__lat       (Structured* This, size_t jlat);
  double atlas__grid__global__Structured__lon       (Structured* This, size_t jlat, size_t jlon);
  void   atlas__grid__global__Structured__lonlat    (Structured* This, size_t jlat, size_t jlon, double crd[]);
  void   atlas__grid__global__Structured__latitudes (Structured* This, const double* &lats, size_t &size);
  int    atlas__grid__global__Structured__reduced   (Structured* This);
}
#undef CONFIG
//------------------------------------------------------------------------------

} // namespace global
} // namespace grid
} // namespace atlas

#endif // Structured_h
