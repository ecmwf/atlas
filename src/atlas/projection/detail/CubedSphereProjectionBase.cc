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
  ATLAS_TRACE( "CubedSphereProjectionBase::CubedSphereProjectionBase" );
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
  ATLAS_TRACE( "CubedSphereProjectionBase::getTile1LonLat" );
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
  ATLAS_TRACE( "CubedSphereProjectionBase::xy2lonlat1" );
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
  double xyz[3];
  double angle;

  this->getTile1LonLat(xytll[0], xytll[1], lonlat);
  ProjectionUtilities::sphericalToCartesian(lonlat, xyz);

  angle = -M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);

  ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

  xytll[3+LON] = lonlat[LON];
  xytll[3+LAT] = lonlat[LAT];
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat3( double xytll[] ) const {
  //  Face 3: rotate -90.0 degrees about z axis
  //          rotate  90.0 degrees about x axis

  double lonlat[2];
  double xyz[3];
  double angle;

  this->getTile1LonLat(xytll[0], xytll[1], lonlat);

  ProjectionUtilities::sphericalToCartesian(lonlat, xyz);

  angle = -M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dX(angle, xyz);

  ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

  xytll[3+LON] = lonlat[LON];
  xytll[3+LAT] = lonlat[LAT];
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat4( double xytll[] ) const {
  //  Face 4: rotate -180.0 degrees about z axis
  //          rotate   90.0 degrees about x axis

  double lonlat[2];
  double xyz[3];
  double angle;

  this->getTile1LonLat(xytll[0], xytll[1], lonlat);

  ProjectionUtilities::sphericalToCartesian(lonlat, xyz);

  angle = -M_PI;
  ProjectionUtilities::rotate3dZ(angle, xyz);
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dX(angle, xyz);

  ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

  xytll[3+LON] = lonlat[LON];
  xytll[3+LAT] = lonlat[LAT];
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat5( double xytll[] ) const {
  //  Face 5: rotate 90.0 degrees about z axis
  //          rotate 90.0 degrees about y axis

  double lonlat[2];
  double xyz[3];
  double angle;

  this->getTile1LonLat(xytll[0], xytll[1], lonlat);

  ProjectionUtilities::sphericalToCartesian(lonlat, xyz);

  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);
  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dY(angle, xyz);

  ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

  xytll[3+LON] = lonlat[LON];
  xytll[3+LAT] = lonlat[LAT];
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::xy2lonlat6( double xytll[] ) const {
  //  Face 6: rotate 90.0 degrees about y axis
  //          rotate  0.0 degrees about y axis

  double lonlat[2];
  double xyz[3];
  double angle;

  this->getTile1LonLat(xytll[0], xytll[1], lonlat);

  ProjectionUtilities::sphericalToCartesian(lonlat, xyz);

  angle = M_PI / 2.0;
  ProjectionUtilities::rotate3dY(angle, xyz);
  angle = 0.0;
  ProjectionUtilities::rotate3dZ(angle, xyz);

  ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

  xytll[3+LON] = lonlat[LON];
  xytll[3+LAT] = lonlat[LAT];
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy1( double xytll[] ) const {
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy2( double xytll[] ) const {
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy3( double xytll[] ) const {
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy4( double xytll[] ) const {
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy5( double xytll[] ) const {
}

// -------------------------------------------------------------------------------------------------

void CubedSphereProjectionBase::latlon2xy6( double xytll[] ) const {
}

// -------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
