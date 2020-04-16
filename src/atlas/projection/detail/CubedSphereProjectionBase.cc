/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereProjectionBase.h"

#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------

CubedSphereProjectionBase::CubedSphereProjectionBase( const eckit::Parametrisation& params )
                                                            : tile1LonsArray_(), tile1LatsArray_() {
  // Get cube sphere face dimension
  params.get("CubeNx", cubeNx_);

  // Arrays to hold projection for (0,0) centered tile
  tile1LonsArray_.reset(new ArrayLatLon_(cubeNx_+1, cubeNx_+1));
  tile1LatsArray_.reset(new ArrayLatLon_(cubeNx_+1, cubeNx_+1));

}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::lonlat2xy( double xytll[] ) const {
  this->lonlat2xyTile.at(xytll[2])(xytll);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat( double xytll[] ) const {
  double xyz[3];
  double lonlat[2];
  xy2lonlatTile.at(xytll[2])(xytll);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::getTile1LonLat(double x, double y, double lonlat[]) const {
  // Get the array view to to tile 1
  auto tile1Lats = getLatArray();
  auto tile1Lons = getLonArray();
  int xi = static_cast<int>(x);
  int yi = static_cast<int>(y);
  lonlat[LON] = tile1Lons(xi, yi);
  lonlat[LAT] = tile1Lats(xi, yi);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat1( double xytll[] ) const {
  //  Face 1, nothing to do.
  double lonlat[2];
  this->getTile1LonLat(xytll[0], xytll[1], lonlat);
  xytll[3+LON] = lonlat[LON];
  xytll[3+LAT] = lonlat[LAT];
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat2( double xytll[] ) const {
  //  Face 2: rotate -90.0 degrees about z axis

  double lonlat[2];
  this->getTile1LonLat(xytll[0], xytll[1], lonlat);

  double xyz[3];
  ProjectionUtilities::sphericalToCartesian(lonlat, xyz);

  double angle = M_PI / 2.0;
  ProjectionUtilities::rotate3d.at(3)(angle, xyz);

  ProjectionUtilities::sphericalToCartesian(xyz, lonlat);

  xytll[3+LON] = lonlat[LON];
  xytll[3+LAT] = lonlat[LAT];
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat3( double xytll[] ) const {
  //  Face 3: rotate -90.0 degrees about z axis
  //          rotate  90.0 degrees about x axis
  std::cout << "xy2lonlat2" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat4( double xytll[] ) const {
  //  Face 4: rotate -180.0 degrees about z axis
  //          rotate   90.0 degrees about x axis
  std::cout << "xy2lonlat3" << std::endl;};

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat5( double xytll[] ) const {
  //  Face 5: rotate 90.0 degrees about z axis
  //          rotate 90.0 degrees about y axis
  std::cout << "xy2lonlat4" << std::endl;};

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat6( double xytll[] ) const {
  //  Face 6: rotate 90.0 degrees about y axis
  //          rotate 90.0 degrees about y axis
  std::cout << "xy2lonlat5" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy1( double xytll[] ) const {
  std::cout << "latlon2xy0" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy2( double xytll[] ) const {
  std::cout << "latlon2xy1" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy3( double xytll[] ) const {
  std::cout << "latlon2xy2" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy4( double xytll[] ) const {
  std::cout << "latlon2xy3" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy5( double xytll[] ) const {
  std::cout << "latlon2xy4" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy6( double xytll[] ) const {
  std::cout << "latlon2xy5" << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
