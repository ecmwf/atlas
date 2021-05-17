/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereEquiAnglProjection.h"

#include <cmath>
#include <iostream>

#include "eckit/config/Parametrisation.h"
#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/projection/detail/ProjectionUtilities.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"


namespace atlas {
namespace projection {
namespace detail {

// -------------------------------------------------------------------------------------------------

CubedSphereEquiAnglProjection::CubedSphereEquiAnglProjection( const eckit::Parametrisation& params )
                                                               : CubedSphereProjectionBase(params) {


 /*
  // Get Base data
  const auto cubeNx = getCubeNx();
  auto tile1Lats = getLatArray();
  auto tile1Lons = getLonArray();

  // Cartesian coordinate starting points
  const double rsq3 = 1.0/sqrt(3.0);

  // Equiangular
  const double dp = 0.5*M_PI/static_cast<double>(cubeNx);

  double xyz[3];
  double lonlat[2];

  for ( int ix = 0; ix < cubeNx+1; ix++ ) {
    for ( int iy = 0; iy < cubeNx+1; iy++ ) {
      // Grid points in cartesian coordinates
      xyz[XX] = -rsq3;
      xyz[YY] = -rsq3*tan(-0.25*M_PI+static_cast<double>(ix)*dp);
      xyz[ZZ] = -rsq3*tan(-0.25*M_PI+static_cast<double>(iy)*dp);

      ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

      if (lonlat[LON] < 0.0) {
        lonlat[LON] += 2.0*M_PI;
      }

      tile1Lons(ix, iy) = lonlat[LON] - M_PI;
      tile1Lats(ix, iy) = lonlat[LAT];
    }
  }
  */
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiAnglProjection::lonlat2xy( double crd[] ) const {

    const double rsq3 = 1.0/sqrt(3.0);
     double xyz[3];
     double ab[2]; // alpha-beta coordinate

   // need lonlat2 alpha beta tile

   //
   // if we can work this out then we can go into Cartesian
   //
   // inverse rotate to tile 0
   //
   // then inverse tan to get alpha beta.

   // convert degrees tp radians
   crd[0] *= M_PI/180.;
   crd[1] *= M_PI/180.;

   // find tile which this lonlat is linked to
   idx_t t = identityTileFromLonLat(crd);

   ProjectionUtilities::sphericalToCartesian(crd, xyz);
   tileRotateInverse.at(t)(xyz);

   //now should be tile 0 - now calculate (alpha, beta) in radians.
   // should be between - pi/4 and pi/4
   ab[0] = - atan2(-xyz[YY], rsq3)/2.0;
   ab[1] = atan2(-xyz[ZZ], rsq3)/2.0;

   std::cout << "lonlat2xy ab : " << ab[0] << " " << ab[1] << std::endl;

   CubedSphereProjectionBase::alphabetatt2xy(t, ab, crd);

  //CubedSphereProjectionBase::lonlat2xy(crd);

}

// -------------------------------------------------------------------------------------------------
// input should be Willems xy coordinate in degrees
//
void CubedSphereEquiAnglProjection::xy2lonlat( double crd[] ) const {

    const double rsq3 = 1.0/sqrt(3.0);
    double xyz[3];
    double ab[2]; // alpha-beta coordinate
    idx_t t;  // tile index
    double lonlat[2];

    // calculate xy (in degrees) to alpha beta (in radians) and t - tile index.
    CubedSphereProjectionBase::xy2alphabetat(crd, t, ab);

    std::cout << "xy2lonlat:: crd t ab  : "  << crd[0] << " " << crd[1] << " " << t << " " << ab[0] << " " << ab[1] << std::endl;

    xyz[0] = -rsq3;
    xyz[1] = -rsq3*tan(ab[0]);
    xyz[2] = -rsq3*tan(ab[1]);

    ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

    if (lonlat[LON] < 0.0) {
      lonlat[LON] += 2.0*M_PI;
    }
    lonlat[LON] = lonlat[LON] - M_PI;

    std::cout << "xy2lonlat:: lonlat before rotation : "  << lonlat[0] << " " << lonlat[1]  << std::endl;

    // Convert to cartesian
    ProjectionUtilities::sphericalToCartesian(lonlat, xyz);

    // Perform tile specific rotation
    tileRotate.at(t)(xyz);

    // Back to latlon
    ProjectionUtilities::cartesianToSpherical(xyz, lonlat);

    // Shift longitude
    /*
    if (shiftLon_ != 0.0) {
      lonlat[LON] = lonlat[LON] + shiftLon_*atlas::util::Constants::degreesToRadians();
      if (lonlat[LON] < -M_PI) {lonlat[LON] =  2*M_PI + lonlat[LON];}
      if (lonlat[LON] >  M_PI) {lonlat[LON] = -2*M_PI + lonlat[LON];}
    }
    */
    // To 0, 360
    if (lonlat[LON] < 0.0) {
      lonlat[LON] = 2*M_PI + lonlat[LON];
    }

    // longitude does not make sense at the poles - set to 0.
    if ( std::abs(std::abs(lonlat[LAT]) - M_PI_2) < 1e-13) lonlat[LON] = 0.;


    crd[LON] = lonlat[LON] * 180.0 / M_PI;
    crd[LAT] = lonlat[LAT] * 180.0 / M_PI;


    std::cout << "end of xy2lonlat" << std::endl;
}


// -------------------------------------------------------------------------------------------------

ProjectionImpl::Jacobian CubedSphereEquiAnglProjection::jacobian(const PointLonLat& ) const {
    ATLAS_NOTIMPLEMENTED;
}


// -------------------------------------------------------------------------------------------------

CubedSphereEquiAnglProjection::Spec CubedSphereEquiAnglProjection::spec() const {
  // Fill projection specification
  Spec proj;
  proj.set( "type", static_type() );
  return proj;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiAnglProjection::hash( eckit::Hash& h ) const {
  // Add to hash
  h.add( static_type() );
  CubedSphereProjectionBase::hash(h);
}

/*
 * function identify_panel(x,y,z) result(panel_id)
524
525	  implicit none
526	  real(kind=r_def), intent(in) :: x,y,z
527	  integer(kind=i_def) :: panel_id
528
529	  real(kind=r_def) :: lon, lat, radius
530
531	  if (z<-abs(x) .and. z<-abs(y))then
532	    panel_id = 6
533	  elseif (z>abs(x) .and. z>abs(y))then
534	    panel_id = 5
535	  else
536	    call xyz2llr(x, y, z, lon, lat, radius)
537	    panel_id = identify_longitude_sector(lon)
538	  end if
539
540	end function identify_panel

202	!------------------------------------------------------------------------------
203	!>  @brief  Converts Cartesian coordinates to longitude, latitude and radius
204	!!  @param[in]   x     Cartesian x coordinate to convert.
205	!!  @param[in]   y     Cartesian y coordinate to convert.
206	!!  @param[in]   z     Cartesian z coordinate to convert.
207	!!  @param[out]  longitude  -PI   <  Longitude <= PI   (radians).
208	!!  @param[out]  latitude   -PI/2 <= Latitude  <= PI/2 (radians).
209	!!  @param[out]  radius     Radius of the sphere(m).
210	!------------------------------------------------------------------------------
211	subroutine xyz2llr( x, y, z, &
212	                    longitude, latitude, radius )
213
214	  implicit none
215
216	  real(r_def), intent(in)  :: x, y, z
217	  real(r_def), intent(out) :: longitude, latitude, radius
218
219	  ! Local variables
220	  real(r_def) :: tan_longitude, tan_latitude
221	  real(r_def) :: tol = 10.0e-8_r_def
222
223
224	  ! Calculate longitude in range
225	  ! -180 < longitude <= 180
226	  if (x == 0.0_r_def) then
227
228	    if (y >= 0.0_r_def) then
229	      longitude =  0.5_r_def*PI
230	    else
231	      longitude = -0.5_r_def*PI
232	    end if
233
234	  else
235
236	    tan_longitude = y/x
237	    longitude = atan(tan_longitude)
238
239	    if (x < 0.0_r_def) then
240	      if (y >= 0.0_r_def) then
241	        longitude = longitude + PI
242	      else
243	        longitude = longitude - PI
244	      end if
245	    end if
246
247	  end if
248
249
250	  ! Calculate latitude in range
251	  ! -90 <= longitude <= +90
252	  radius = sqrt(x*x + y*y)
253	  if ( abs(radius) < tol ) then
254	    if (z > 0.0_r_def ) then
255	      latitude =  0.5_r_def*PI
256	    else
257	      latitude = -0.5_r_def*PI
258	    end if
259	    ! Ensure consisent value for longitude is
260	    ! output for Latitudes of -90/90.
261	    longitude = 0.0_r_def
262	  else
263	    tan_latitude = z/radius
264	    latitude = atan(tan_latitude)
265	  end if
266
267	  ! Radius
268	  radius = sqrt(x*x + y*y + z*z)
269
270	end subroutine xyz2llr

*/



// assuming that crd[] here holds the latitude/longitude in RADIANS.
// expecting longitude range between 0 and 2* PI
idx_t CubedSphereEquiAnglProjection::identityTileFromLonLat(const double crd[]) const {
    idx_t t(-1);
    double xyz[3];

    ProjectionUtilities::sphericalToCartesian(crd, xyz, false, true);

    const double cornerLat= std::asin(1./std::sqrt(3.0));
      // the magnitude of the latitude at the corners of the cube (not including the sign)
      // in radians.

    const double & lon = crd[LON];
    const double & lat = crd[LAT];

    double zPlusAbsX = xyz[2] + abs(xyz[0]);
    double zPlusAbsY = xyz[2] + abs(xyz[1]);
    double zMinusAbsX = xyz[2] - abs(xyz[0]);
    double zMinusAbsY = xyz[2] - abs(xyz[1]);

    std::cout << "lon lat radians " << lon << " " << lat  << " " << abs(lat - cornerLat) << std::endl;

    if (lon >= 1.75 * M_PI  || lon < 0.25 * M_PI) {
        if  ( (zPlusAbsX <= 0.) && (zPlusAbsY <= 0.) ) {
           t = 2;
        } else if ( (zMinusAbsX > 0.) && (zMinusAbsY > 0.) ) {
           t = 5;
        } else {
           t = 0;
        }
        // extra point corner point
        std::cout << "first extra corner " << abs(lon - 1.75 * M_PI) << " " << abs(lat - cornerLat) << std::endl;
        if ( (abs(lon - 1.75 * M_PI) < 1e-13) &&
             (abs(lat - cornerLat) < 1e-13) ) t = 0;
    }

    if (lon >= 0.25 * M_PI  && lon < 0.75 * M_PI) {
        // interior
        if  ( (zPlusAbsX <= 0.) && (zPlusAbsY <= 0.) ) {
            t = 2;
        } else if  ( (zMinusAbsX > 0.) && (zMinusAbsY > 0.) ) {
            t = 5;
        } else {
            t = 1;
        }

    }

    if (lon >= 0.75 * M_PI  && lon < 1.25 * M_PI) {
        // interior
        if  ( (zPlusAbsX < 0.) && (zPlusAbsY < 0.) ) {
            t = 2;
        } else if  ( (zMinusAbsX >= 0.) && (zMinusAbsY >= 0.) ) {
            t = 5;
        } else {
            t = 3;
        }
        // extra point corner point
        std::cout << "second extra corner " << abs(lon - 0.75 * M_PI) << " " << abs(lat + cornerLat) << std::endl;
        if ( (abs(lon - 0.75 * M_PI) < 1e-13) &&
             (abs(lat + cornerLat) < 1e-13) ) t = 1;
    }

    if (lon >= 1.25 * M_PI  && lon < 1.75 * M_PI) {
        // interior
        if  ( (zPlusAbsX < 0.) && (zPlusAbsY <=0.) ) {
            t = 2;
        } else if  ( (zMinusAbsX >= 0.) && (zMinusAbsY >= 0.) ) {
            t = 5;
        } else {
            t = 4;
        }
    }

    std::cout << "identifyTileFromLonLat:: lonlat  ...t " <<
        crd[0] << " " << crd[1]  << " " <<  0.25 * M_PI  << " " << 1.75 * M_PI << " " << t << std::endl;
    std::cout << "identifyTileFromLonLat::  xyz, zPlusAbsX ... " <<
        xyz[0] << " " << xyz[1]  << " " << xyz[2] << " " <<
        zPlusAbsX << " " << zPlusAbsY << " " <<
        zMinusAbsX << " " << zMinusAbsY << " " <<  std::endl;


    /*

        // call xyz2llr(x, y, z, lon, lat, radius)
        const double & lon = crd[LON];

        if (lon >= -M_PI_4 && lon< M_PI_4) {
           t = 0;
        } else if (lon >= M_PI_4 && lon< 3.0 * M_PI_4) {
           t = 1;


        if ((lon > -5.0 * M_PI_4) && lon < -3.0 * M_PI_4) {
            t = 3;
        } else if (lon > -3.0 * M_PI_4 && lon < -M_PI_4) {
            t = 4;
        } else if (lon > -M_PI_4 && lon< M_PI_4) {
            t = 0;
        } else if (lon > M_PI_4 && lon< 3.0 * M_PI_4) {
            t = 1;
        } else if (lon > 3.0 * M_PI_4 && lon < 5.0 * M_PI_4) {
            t = 3;
        } else if (lon >  5.0 * M_PI_4 && lon < 7.0 * M_PI_4) {
            t = 4;
        } else if (lon > 7.0 * M_PI_4) {
            t = 0;
        }
    }
    // edges


    // corners
   */

    return t;
}

// -------------------------------------------------------------------------------------------------

namespace {
static ProjectionBuilder<CubedSphereEquiAnglProjection>
       register_1( CubedSphereEquiAnglProjection::static_type() );
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
