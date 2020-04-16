/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <functional>

#include "eckit/config/Parametrisation.h"

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

    // projection and inverse projection
    void xy2lonlat( double crd[] ) const;
    void lonlat2xy( double crd[] ) const;

    // Get the latlon for tile 1
    void getTile1LonLat(double, double, double[]) const;

    // Functions for xy to latlon on each tile
    void xy2lonlat1( double[] ) const;
    void xy2lonlat2( double[] ) const;
    void xy2lonlat3( double[] ) const;
    void xy2lonlat4( double[] ) const;
    void xy2lonlat5( double[] ) const;
    void xy2lonlat6( double[] ) const;

    std::vector<std::function<void(double[])>> xy2lonlatTile =
      {[this](double xytll[]){this->xy2lonlat1(xytll);},
       [this](double xytll[]){this->xy2lonlat2(xytll);},
       [this](double xytll[]){this->xy2lonlat3(xytll);},
       [this](double xytll[]){this->xy2lonlat4(xytll);},
       [this](double xytll[]){this->xy2lonlat5(xytll);},
       [this](double xytll[]){this->xy2lonlat6(xytll);}
      };

    // Functions for latlon to xy on each tile
    void latlon2xy1( double[] ) const;
    void latlon2xy2( double[] ) const;
    void latlon2xy3( double[] ) const;
    void latlon2xy4( double[] ) const;
    void latlon2xy5( double[] ) const;
    void latlon2xy6( double[] ) const;

    std::vector<std::function<void(double[])>> lonlat2xyTile =
      {[this](double xytll[]){this->latlon2xy1(xytll);},
       [this](double xytll[]){this->latlon2xy2(xytll);},
       [this](double xytll[]){this->latlon2xy3(xytll);},
       [this](double xytll[]){this->latlon2xy4(xytll);},
       [this](double xytll[]){this->latlon2xy5(xytll);},
       [this](double xytll[]){this->latlon2xy6(xytll);}
      };

    // Array views for accessing data of tile 0 projection
    ArrayViewLatLon_ getLatArray() const {
      ATLAS_TRACE( "CubedSphereProjectionBase::getLatArray" );
      array::ArrayView<double, 2> tile1Lats = array::make_view<double, 2>( *tile1LatsArray_ );
      return tile1Lats;
    };
    ArrayViewLatLon_ getLonArray() const {
      ATLAS_TRACE( "CubedSphereProjectionBase::getLonArray" );
      array::ArrayView<double, 2> tile1Lats = array::make_view<double, 2>( *tile1LonsArray_ );
      return tile1Lats;
    };

    int getCubeNx() const { return cubeNx_; };

  private:
    int cubeNx_;
    std::unique_ptr<ArrayLatLon_> tile1LonsArray_;
    std::unique_ptr<ArrayLatLon_> tile1LatsArray_;
};

}  // namespace detail
}  // namespace projection
}  // namespace atlas
