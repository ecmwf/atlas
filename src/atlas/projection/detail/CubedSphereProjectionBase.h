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

  protected:
    idx_t tile( const double xy[] ) const {
      // Assume one face-edge is of length 90 degrees.
      //
      //   y ^
      //     |
      //    135           -------
      //     |           |       |
      //     |           |   2   |
      //     |           |       |
      //     45   ------- ------- ------- -------
      //     |   |       |       |       |       |
      //     |   |   0   |   1   |   3   |   4   |
      //     |   |       |       |       |       |
      //    -45   -------  ------ ------- -------
      //     |                           |       |
      //     |                           |   5   |
      //     |                           |       |
      //   -135                           -------
      //     ----0-------90------180-----270----360--->  x
      const double x = xy[0]/90.;
      const double& y = xy[1];
      if( x < 2. ) {
        return  y > 45. ? 2 : std::floor(x);
      }
      else {
        return y < -45. ? 5 : std::floor(x+1.);
      }
  }

    void xy2xyt(const double xy[], double xyt[]) const {
        // xy is in degrees while xyt is in radians
        // (alpha, beta) and tiles.
        double normalisedX = xy[0]/90.;
        double normalisedY = (xy[1] + 135.)/90.;
        xyt[0] = (normalisedX - std::floor(normalisedX))* M_PI_2 - M_PI_4;
        xyt[1] = (normalisedY - std::floor(normalisedY))* M_PI_2 - M_PI_4;
        xyt[2] = tile(xy);
    }

    void xyt2xy(const double xyt[], double xy[]) const {
        // xy is in degrees while xyt is in radians
        // (alpha, beta) and tiles.
        std::vector<double> xOffset{0., 90., 90., 180, 270, 270};
        std::vector<double> yOffset{-45., -45, 45, -45, -45, -135};
        double normalisedX = (xyt[0] + M_PI_4)/M_PI_2;
        double normalisedY = (xyt[1] + M_PI_4)/M_PI_2;
        xy[0] = normalisedX * 90. + xOffset[xyt[2]];
        xy[1] = normalisedY * 90. + yOffset[xyt[2]];
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
