/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

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
    void xy2lonlat( double xytll[] ) const;
    void lonlat2xy( double llxyt[] ) const;

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
      array::ArrayView<double, 2> tile1Lats = array::make_view<double, 2>( *tile1LonsArray_ );
      return tile1Lats;
    }

    int getCubeNx() const { return cubeNx_; }

    void schmidtTransform(double, double, double, double[]) const;

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
