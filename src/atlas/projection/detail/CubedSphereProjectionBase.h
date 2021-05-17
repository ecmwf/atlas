/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <cmath>
#include <functional>

#include "eckit/config/Parametrisation.h"
#include "eckit/utils/Hash.h"

#include "atlas/array.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"

namespace atlas {
namespace projection {
namespace detail {

class CubedSphereProjectionBase {
  typedef array::ArrayT<double> ArrayLatLon_;
  typedef array::ArrayView<double, 2> ArrayViewLatLon_;
  public:
    // constructor
    CubedSphereProjectionBase( const eckit::Parametrisation& );

    void hash( eckit::Hash& ) const;

    // projection and inverse projection
    void xy2lonlat( double crd[] ) const;
    void lonlat2xy( double crd[] ) const;

    // Functions for xy to latlon on each tile
    void tile1Rotate( double[] ) const;
    void tile2Rotate( double[] ) const;
    void tile3Rotate( double[] ) const;
    void tile4Rotate( double[] ) const;
    void tile5Rotate( double[] ) const;
    void tile6Rotate( double[] ) const;

    std::vector<std::function<void(double[])>> tileRotate =
      {[this](double xyz[]){this->tile1Rotate(xyz);},
       [this](double xyz[]){this->tile2Rotate(xyz);},
       [this](double xyz[]){this->tile3Rotate(xyz);},
       [this](double xyz[]){this->tile4Rotate(xyz);},
       [this](double xyz[]){this->tile5Rotate(xyz);},
       [this](double xyz[]){this->tile6Rotate(xyz);}
      };

    // Functions for latlon to xy on each tile
    void tile1RotateInverse( double[] ) const;
    void tile2RotateInverse( double[] ) const;
    void tile3RotateInverse( double[] ) const;
    void tile4RotateInverse( double[] ) const;
    void tile5RotateInverse( double[] ) const;
    void tile6RotateInverse( double[] ) const;

    std::vector<std::function<void(double[])>> tileRotateInverse =
      {[this](double xyz[]){this->tile1RotateInverse(xyz);},
       [this](double xyz[]){this->tile2RotateInverse(xyz);},
       [this](double xyz[]){this->tile3RotateInverse(xyz);},
       [this](double xyz[]){this->tile4RotateInverse(xyz);},
       [this](double xyz[]){this->tile5RotateInverse(xyz);},
       [this](double xyz[]){this->tile6RotateInverse(xyz);}
      };

    // Array views for accessing data of tile 0 projection
    ArrayViewLatLon_ getLatArray() const {
      ATLAS_TRACE( "CubedSphereProjectionBase::getLatArray" );
      array::ArrayView<double, 2> tile1Lats = array::make_view<double, 2>( *tile1LatsArray_ );
      return tile1Lats;
    }
    ArrayViewLatLon_ getLonArray() const {
      ATLAS_TRACE( "CubedSphereProjectionBase::getLonArray" );
      array::ArrayView<double, 2> tile1Lons = array::make_view<double, 2>( *tile1LonsArray_ );
      return tile1Lons;
    }

    int getCubeNx() const { return cubeNx_; }

    void schmidtTransform(double, double, double, double[]) const;

    void xy2alphabetat(const double xy[], idx_t & t, double ab[]) const {
        // xy is in degrees while ab is in radians
        // ab are the  (alpha, beta) coordinates and t is the tile index.
        t = tile(xy);
        std::vector<double> xOffset{0., 1., 1., 2., 3., 3.};
        std::vector<double> yOffset{1., 1., 2., 1., 1., 0.};

        double normalisedX = xy[0]/90.;
        double normalisedY = (xy[1] + 135.)/90.;
        ab[0] = (normalisedX - xOffset[t])* M_PI_2 - M_PI_4;
        ab[1] = (normalisedY - yOffset[t])* M_PI_2 - M_PI_4;

    }

    void alphabetatt2xy(const idx_t & t, const double ab[],  double xy[]) const {
        // xy is in degrees while ab is in radians
        // (alpha, beta) and tiles.
        std::vector<double> xOffset{0., 90., 90., 180, 270, 270};
        std::vector<double> yOffset{-45., -45, 45, -45, -45, -135};
        double normalisedX = (ab[0] + M_PI_4)/M_PI_2;
        double normalisedY = (ab[1] + M_PI_4)/M_PI_2;
        xy[0] = normalisedX * 90. + xOffset[t];
        xy[1] = normalisedY * 90. + yOffset[t];
  }

  protected:
    idx_t tile( const double xy[] ) const {
      // Assume one face-edge is of length 90 degrees.
      //
      //   y ^
      //     |
      //    135              ----------
      //     |              |     ^    |
      //     |              |          |
      //     |              |=<   2   <|
      //     |              |     v    |
      //     |              |     =    |
      //     45  0----------2----------3----------4----------
      //     |   |    ^     |     ^    |    =     |     =    |
      //     |   |          |          |    ^     |     ^    |
      //     |   |=<  0    <|=<   1   <|=<  3    <|=<   4   <|
      //     |   |    v     |     v    |          |          |
      //     |   |    =     |     =    |    v     |     v    |
      //    -45  0 ---------1----------1----------5----------
      //     |                                    |     =    |
      //     |                                    |     ^    |
      //     |                                    |=<   5   <|
      //     |                                    |          |
      //     |                                    |     v    |
      //   -135                                    ----------(5 for end iterator)
      //     ----0---------90--------180--------270--------360--->  x
      /*
      const double x = xy[0]/90.;
      const double& y = xy[1];
      if( x < 2. ) {
        return  y > 45. ? 2 : std::floor(x);
      }
      else {
        return y < -45. ? 5 : std::floor(x+1.);
      }
      */
      idx_t t{-1};

      if ((xy[0] >= 0.) && ( xy[1] >= -45.) && (xy[0] < 90.) && (xy[1] < 45.)) {
         t = 0;
      } else if ((xy[0] >= 90.) && ( xy[1] >= -45.) && (xy[0] < 180.) && (xy[1] < 45.)) {
         t = 1;
      } else if ((xy[0] >= 90.) && ( xy[1] >= 45.) && (xy[0] < 180.) && (xy[1] < 135.)) {
         t = 2;
      } else if ((xy[0] >= 180.) && ( xy[1] > -45.) && (xy[0] < 270.) && (xy[1] <= 45.)) {
         t = 3;
      } else if ((xy[0] >= 270.) && ( xy[1] > -45.) && (xy[0] < 360.) && (xy[1] <= 45.)) {
         t = 4;
      } else if ((xy[0] >= 270.) && ( xy[1] > -135.) && (xy[0] < 360.) && (xy[1] <= -45.)) {
         t = 5;
      }

      // extra points
      if ((xy[0] == 0.) && (xy[1] == 45.)) t = 0;
      if ((xy[0] == 180.) && (xy[1] == -45.)) t = 1;

      // for end iterator !!!!
      if ((xy[0] == 360.) && (xy[1] == -135.)) t = 5;

      return t;

  }

  private:
    int cubeNx_;
    std::unique_ptr<ArrayLatLon_> tile1LonsArray_;
    std::unique_ptr<ArrayLatLon_> tile1LatsArray_;
    // Shift entire grid
    double shiftLon_;
    // Schmidt transform
    bool doSchmidt_;
    double stretchFac_;
    double targetLon_;
    double targetLat_;
};

}  // namespace detail
}  // namespace projection
}  // namespace atlas
