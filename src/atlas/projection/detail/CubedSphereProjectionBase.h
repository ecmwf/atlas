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
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace projection {
namespace detail {

class CubedSphereProjectionBase : public ProjectionImpl {
  typedef array::ArrayT<double> ArrayLatLon_;
  typedef array::ArrayView<double, 2> ArrayViewLatLon_;
  public:
    // constructor
    CubedSphereProjectionBase( const eckit::Parametrisation& );

    void hash( eckit::Hash& ) const;

    // projection and inverse projection
    void xy2lonlatpost( double xyz[], const idx_t & t, double crd[] ) const;
    void lonlat2xypre( double crd[], idx_t & t, double xyz[] ) const;

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

    void schmidtTransform(double, double, double, double[]) const;

    void xy2alphabetat(const double xy[], idx_t & t, double ab[]) const {
        // xy is in degrees while ab is in radians
        // ab are the  (alpha, beta) coordinates and t is the tile index.
        t = tileFromXY(xy);
        std::vector<double> xOffset{0., 1., 1., 2., 3., 3.};
        std::vector<double> yOffset{1., 1., 2., 1., 1., 0.};

        double normalisedX = xy[XX]/90.;
        double normalisedY = (xy[YY] + 135.)/90.;
        ab[LON] = (normalisedX - xOffset[t])* M_PI_2 - M_PI_4;
        ab[LAT] = (normalisedY - yOffset[t])* M_PI_2 - M_PI_4;

    }

    void alphabetatt2xy(const idx_t & t, const double ab[],  double xy[]) const {
        // xy is in degrees while ab is in radians
        // (alpha, beta) and tiles.
        std::vector<double> xOffset{0., 90., 90., 180, 270, 270};
        std::vector<double> yOffset{-45., -45, 45, -45, -45, -135};
        double normalisedX = (ab[LON] + M_PI_4)/M_PI_2;
        double normalisedY = (ab[LAT] + M_PI_4)/M_PI_2;
        xy[XX] = normalisedX * 90. + xOffset[t];
        xy[YY] = normalisedY * 90. + yOffset[t];
    }

    idx_t tileFromXY(const double xy[]) const;

  protected:  

    idx_t tileFromLonLat(const double crd[]) const;

  private:
    // Shift entire grid
    double shiftLon_;
    // Schmidt transform
    bool doSchmidt_;
    double stretchFac_;
    double targetLon_;
    double targetLat_;
    const double epsilon_;
};

}  // namespace detail
}  // namespace projection
}  // namespace atlas
