/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereEquiDistFV3Projection.h"

#include <cmath>
#include <iostream>
#include <iomanip>

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

CubedSphereEquiDistFV3Projection::CubedSphereEquiDistFV3Projection( const eckit::Parametrisation& params )
                                                                  : CubedSphereProjectionBase(params) {

  //  Equidistant cubedsphere used in the FV3 dynamical core. Equidistant along the 4 edges of the cubed sphere
  //  Properties:
  //   - defined by intersections of great circles
  //   - max(dx,dy; global) / min(dx,dy; global) = sqrt(2) = 1.4142
  //   - Max(aspect ratio) = 1.06089
  //   - the N-S coordinate curves are const longitude on the 4 faces with equator


  // Get Base data
  const auto cubeNx = getCubeNx();
  auto tile1Lats = getLatArray();
  auto tile1Lons = getLonArray();

  // Indices
  int ix;
  int iy;

  // Array for cartesian coordinates
  array::ArrayT<double> xyzArrayT( 3, cubeNx+1, cubeNx+1 );
  array::ArrayView<double, 3> xyzArray = array::make_view<double, 3>( xyzArrayT );

  // Zero the arrays
  for (ix = 0; ix < cubeNx + 1; ix++ ) {
    for (iy = 0; iy < cubeNx + 1; iy++ ) {
      tile1Lons(ix, iy) = 0.0;
      tile1Lats(ix, iy) = 0.0;
    }
  }
  for (int n = 0; n < 3; n++ ) {
    for (ix = 0; ix < cubeNx + 1; ix++ ) {
      for (iy = 0; iy < cubeNx + 1; iy++ ) {
        xyzArray(n, ix, iy) = 0.0;
      }
    }
  }

  // Cartesian coordinate starting points
  const double rsq3 = 1.0/sqrt(3.0);
  const double alpha = asin( rsq3 );
  const double dely = 2.0*alpha / static_cast<double>(cubeNx);

  // Define east and west edges
  for ( iy = 0; iy < cubeNx+1; iy++ ) {
    tile1Lons(0,      iy) = 0.75*M_PI;                              // West edge
    tile1Lons(cubeNx, iy) = 1.25*M_PI;                              // East edge
    tile1Lats(0,      iy) = -alpha + dely*static_cast<double>(iy);  // West edge
    tile1Lats(cubeNx, iy) = tile1Lats(0, iy);                       // East edge
  }

  // Get North-South edges by symmetry
  double lonlat1[2];
  double lonlat2[2];
  double lonlat3[2];
  double lonlat4[2];
  for ( ix = 1; ix < cubeNx; ix++ ) {
    lonlat1[LON] = tile1Lons(0, 0);
    lonlat1[LAT] = tile1Lats(0, 0);
    lonlat2[LON] = tile1Lons(cubeNx, cubeNx);
    lonlat2[LAT] = tile1Lats(cubeNx, cubeNx);
    lonlat3[LON] = tile1Lons(0, ix);
    lonlat3[LAT] = tile1Lats(0, ix);
    ProjectionUtilities::mirror_latlon(lonlat1, lonlat2, lonlat3, lonlat4);
    tile1Lons(ix, 0) = lonlat4[LON];
    tile1Lats(ix, 0) = lonlat4[LAT];
    tile1Lons(ix, cubeNx) = lonlat4[LON];
    tile1Lats(ix, cubeNx) = -lonlat4[LAT];
  }

  // Set 4 corners
  double lonlat_tmp[2];
  double xyz_tmp[3];

  int index_0[4] = {0, cubeNx, 0, cubeNx};
  int index_1[4] = {0, 0, cubeNx, cubeNx};

  for (int n = 0; n < 4; n++) {
    lonlat_tmp[LON] = tile1Lons(index_0[n], index_1[n]);
    lonlat_tmp[LAT] = tile1Lats(index_0[n], index_1[n]);
    ProjectionUtilities::sphericalToCartesian(lonlat_tmp, xyz_tmp, true, true);
    xyzArray(XX, index_0[n], index_1[n]) = xyz_tmp[XX];
    xyzArray(YY, index_0[n], index_1[n]) = xyz_tmp[YY];
    xyzArray(ZZ, index_0[n], index_1[n]) = xyz_tmp[ZZ];
  }


  // Map edges on the sphere back to cube. Intersections at x=-rsq3
  ix = 0;
  for (iy = 1; iy < cubeNx; iy++) {
    lonlat_tmp[LON] = tile1Lons(ix, iy);
    lonlat_tmp[LAT] = tile1Lats(ix, iy);
    ProjectionUtilities::sphericalToCartesian(lonlat_tmp, xyz_tmp, true, true);
    xyzArray(XX, ix, iy) = xyz_tmp[XX];
    xyzArray(YY, ix, iy) = xyz_tmp[YY];
    xyzArray(ZZ, ix, iy) = xyz_tmp[ZZ];

    xyzArray(YY, ix, iy) = -xyzArray(YY, ix, iy)*rsq3/xyzArray(XX, ix, iy);
    xyzArray(ZZ, ix, iy) = -xyzArray(ZZ, ix, iy)*rsq3/xyzArray(XX, ix, iy);
  }

  iy = 0;
  for (ix = 1; ix < cubeNx; ix++) {
    lonlat_tmp[LON] = tile1Lons(ix, iy);
    lonlat_tmp[LAT] = tile1Lats(ix, iy);
    ProjectionUtilities::sphericalToCartesian(lonlat_tmp, xyz_tmp, true, true);
    xyzArray(XX, ix, iy) = xyz_tmp[XX];
    xyzArray(YY, ix, iy) = xyz_tmp[YY];
    xyzArray(ZZ, ix, iy) = xyz_tmp[ZZ];

    xyzArray(YY, ix, iy) = -xyzArray(YY, ix, iy)*rsq3/xyzArray(XX, ix, iy);
    xyzArray(ZZ, ix, iy) = -xyzArray(ZZ, ix, iy)*rsq3/xyzArray(XX, ix, iy);
  }

  for (ix = 0; ix < cubeNx+1; ix++) {
    for (iy = 0; iy < cubeNx+1; iy++) {
      xyzArray(XX, ix, iy) = -rsq3;
    }
  }

  for (ix = 1; ix < cubeNx+1; ix++) {
    for (iy = 1; iy < cubeNx+1; iy++) {
      // Copy y-z face of the cube along iy=1
      xyzArray(YY,ix,iy) = xyzArray(YY,ix,0);
      // Copy along ix=1
      xyzArray(ZZ,ix,iy) = xyzArray(ZZ,0,iy);
    }
  }

  // To lonlat
  for (ix = 0; ix < cubeNx+1; ix++) {
    for (iy = 0; iy < cubeNx+1; iy++) {
      xyz_tmp[XX] = xyzArray(XX, ix, iy);
      xyz_tmp[YY] = xyzArray(YY, ix, iy);
      xyz_tmp[ZZ] = -xyzArray(ZZ, ix, iy);
      ProjectionUtilities::cartesianToSpherical(xyz_tmp, lonlat_tmp);
      if (lonlat_tmp[LON] < 0.0) {
        lonlat_tmp[LON] += 2.0*M_PI;
      }
      tile1Lons(ix, iy) = lonlat_tmp[LON];
      tile1Lats(ix, iy) = lonlat_tmp[LAT];
    }
  }


  // Symmetry
  // --------

  // Fix internal lons to same value at all lats
  for (ix = 1; ix < cubeNx; ix++) {
    for (iy = 1; iy < cubeNx+1; iy++) {
      tile1Lons(ix, iy) = tile1Lons(ix, 0);
    }
  }

  int ip;
  double avg;
  for (ix = 0; ix < cubeNx/2; ix++) {
    for (iy = 0; iy < cubeNx+1; iy++) {
      ip = cubeNx - ix;
      avg = 0.5*(tile1Lons(ix, iy)-tile1Lons(ip, iy));
      tile1Lons(ix, iy) = avg + M_PI;
      tile1Lons(ip, iy) = M_PI - avg;
      avg = 0.5*(tile1Lats(ix, iy)+tile1Lats(ip, iy));
      tile1Lats(ix, iy) = avg;
      tile1Lats(ip, iy) = avg;
    }
  }

  int jp;
  for (ix = 1; ix < cubeNx; ix++) {
    for (iy = 0; iy < cubeNx/2; iy++) {
      jp = cubeNx - iy;
      avg = 0.5*(tile1Lons(ix, iy)+tile1Lons(ix, jp));
      tile1Lons(ix, iy) = avg;
      tile1Lons(ix, jp) = avg;
      avg = 0.5*(tile1Lats(ix, iy)-tile1Lats(ix, jp));
      tile1Lats(ix, iy) = avg;
      tile1Lats(ix, jp) = -avg;
    }
  }

  // Minus pi
  for (ix = 0; ix < cubeNx+1; ix++) {
    for (iy = 0; iy < cubeNx+1; iy++) {
      tile1Lons(ix, iy) = tile1Lons(ix, iy) - M_PI;
    }
  }

  // Force consistency
  double x1;
  double y1;
  int cubeHalf = ceil((static_cast<double>(cubeNx)+1.0)/2.0);

  for (ix = 0; ix < cubeHalf + 1; ix++ ) {
    for (iy = 0; iy < cubeHalf + 1; iy++ ) {

      x1 = 0.25 * ( abs(tile1Lons(ix       , iy       )) + abs(tile1Lons(cubeNx-ix, iy       )) +
                    abs(tile1Lons(ix       , cubeNx-iy)) + abs(tile1Lons(cubeNx-ix, cubeNx-iy)) );

      tile1Lons(ix       , iy       ) = copysign(x1, tile1Lons(ix       ,iy       ));
      tile1Lons(cubeNx-ix, iy       ) = copysign(x1, tile1Lons(cubeNx-ix,iy       ));
      tile1Lons(ix       , cubeNx-iy) = copysign(x1, tile1Lons(ix       ,cubeNx-iy));
      tile1Lons(cubeNx-ix, cubeNx-iy) = copysign(x1, tile1Lons(cubeNx-ix,cubeNx-iy));

      y1 = 0.25 * ( abs(tile1Lats(ix       , iy       )) + abs(tile1Lats(cubeNx-ix, iy       )) +
                    abs(tile1Lats(ix       , cubeNx-iy)) + abs(tile1Lats(cubeNx-ix, cubeNx-iy)) );

      tile1Lats(ix       , iy       ) = copysign(y1, tile1Lats(ix       , iy       ));
      tile1Lats(cubeNx-ix, iy       ) = copysign(y1, tile1Lats(cubeNx-ix, iy       ));
      tile1Lats(ix       , cubeNx-iy) = copysign(y1, tile1Lats(ix       , cubeNx-iy));
      tile1Lats(cubeNx-ix, cubeNx-iy) = copysign(y1, tile1Lats(cubeNx-ix, cubeNx-iy));

      //  Force dateline/greenwich-meridion consitency
      if (cubeNx % 2 != 0) {
        if ( (ix==1+((cubeNx+1)-1)/2.0) ) {
          tile1Lons(ix, iy        ) = 0.0;
          tile1Lons(ix, (cubeNx+1)-(iy-1)) = 0.0;
        }
      }

    }
  }

}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistFV3Projection::lonlat2xy( double crd[] ) const {
  CubedSphereProjectionBase::lonlat2xy(crd);
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistFV3Projection::xy2lonlat( double crd[] ) const {
  CubedSphereProjectionBase::xy2lonlat(crd);
}

// -------------------------------------------------------------------------------------------------

CubedSphereEquiDistFV3Projection::Spec CubedSphereEquiDistFV3Projection::spec() const {
  // Fill projection specification
  Spec proj;
  proj.set( "type", static_type() );
  return proj;
}

// -------------------------------------------------------------------------------------------------

void CubedSphereEquiDistFV3Projection::hash( eckit::Hash& h ) const {
  // Add to hash
  h.add( static_type() );
}

// -------------------------------------------------------------------------------------------------

namespace {
static ProjectionBuilder<CubedSphereEquiDistFV3Projection>
       register_1( CubedSphereEquiDistFV3Projection::static_type() );
}

}  // namespace detail
}  // namespace projection
}  // namespace atlas
