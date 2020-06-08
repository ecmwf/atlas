/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphere.h"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <numeric>

#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"
#include "eckit/utils/Translator.h"

#include "atlas/domain/Domain.h"
#include "atlas/grid/CubedSphereGrid.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/grid/detail/grid/GridFactory.h"
#include "atlas/grid/detail/spacing/CustomSpacing.h"
#include "atlas/grid/detail/spacing/LinearSpacing.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Point.h"
#include "atlas/util/UnitSphere.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {


static eckit::Translator<std::string, int> to_int;

static Domain domain( const Grid::Config& grid ) {
  Grid::Config config;
  if ( grid.get( "domain", config ) ) {
      return Domain( config );
  }
  return Domain();
}

std::string CubedSphere::static_type() {
  return "cubedsphere";
}

std::string CubedSphere::name() const {
  return name_;
}

CubedSphere::CubedSphere( int N, Projection p ) : CubedSphere( CubedSphere::static_type(), N, p ) {}

CubedSphere::CubedSphere( const std::string& name, int CubeNx, Projection projection ) :
                           Grid(), CubeNx_(CubeNx), name_( name )  { // Number of tiles hardwired to 6 at the moment. Regional may need 1
  // Copy members
  projection_ = projection ? projection : Projection();

  // Domain
  domain_ = computeDomain();

  // x and y are the position in a 2D plane for the unfolded cubed-sphere grid, shown in the
  // comments in grid/CubedSphereGrid.h. In order to locate the position in this xy array the start
  // position for each face (tile) of the cube is needed. xs represents the x start position and ys
  // the y start position. Tile 3, 4 and 5 are rotated and ysr provides the start point for y after
  // these rotations.

  xs_[0] = 0*CubeNx;
  xs_[1] = 1*CubeNx;
  xs_[2] = 1*CubeNx;
  xs_[3] = 2*CubeNx;
  xs_[4] = 3*CubeNx;
  xs_[5] = 3*CubeNx;

  ys_[0] = 1*CubeNx;
  ys_[1] = 1*CubeNx;
  ys_[2] = 2*CubeNx;
  ys_[3] = 1*CubeNx+1;
  ys_[4] = 1*CubeNx+1;
  ys_[5] = 0*CubeNx+1;

  ysr_[0] = ys_[0];
  ysr_[1] = ys_[1];
  ysr_[2] = ys_[2];
  ysr_[3] = 2*CubeNx;
  ysr_[4] = 2*CubeNx;
  ysr_[5] = 1*CubeNx;

  // Number of grid points on each face of the tile.
  npts_.push_back(CubeNx*CubeNx+1); // An extra point lies on tile 1
  npts_.push_back(CubeNx*CubeNx+1); // An extra point lies on tile 2
  npts_.push_back(CubeNx*CubeNx);
  npts_.push_back(CubeNx*CubeNx);
  npts_.push_back(CubeNx*CubeNx);
  npts_.push_back(CubeNx*CubeNx);
}

// Provide the domain for the cubed-sphere grid, which is global.
Domain CubedSphere::computeDomain() const {return GlobalDomain(); }

// Destructor
CubedSphere::~CubedSphere() = default;

// Print the name of the Grid
void CubedSphere::print( std::ostream& os ) const {
  os << "CubedSphere(Name:" << name() << ")";
}

// Return the type of this Grid
std::string CubedSphere::type() const {
  return static_type();
}

// Provide a unique identification hash for the grid and the projection.
void CubedSphere::hash( eckit::Hash& h ) const {
  h.add("CubedSphere");
  h.add(CubeNx_);

  // also add projection information
  projection().hash( h );

  // also add domain information, even though already encoded in grid.
  domain().hash( h );
}

// Return the bounding box for the grid, global
RectangularLonLatDomain CubedSphere::lonlatBoundingBox() const {
  return projection_ ? projection_.lonlatBoundingBox( computeDomain() ) : domain();
}

// Return the specification for the grid.
Grid::Spec CubedSphere::spec() const {
  Grid::Spec grid_spec;

  if ( name() == "cubedsphere" ) {
    grid_spec.set( "type", type() );
  }
  else {
    grid_spec.set( "name", name() );
  }
  grid_spec.set( "projection", projection().spec() );
  return grid_spec;
}

// --------------------------------------------------------------------

namespace {
  GridFactoryBuilder<CubedSphere> __register_CubedSphere( CubedSphere::static_type() );
}

// -------------------------------------------------------------------------------------------------

// Specialization based on type of projection
// ------------------------------------------

static class cubedsphere_equiangular : public GridBuilder {
public:
  cubedsphere_equiangular() : GridBuilder( "cubedsphere_equiangular", {"^[Cc][Ss][_-][Ee][Aa][-_]([0-9]+)$"}, {"CSEA<cubedsphere>"} ) {}

  void print( std::ostream& os ) const override {
    os << std::left << std::setw( 20 ) << "CS-EA-<FaceNx>" << "Cubed sphere, equiangular";
  }

  // Factory constructor
  const atlas::Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
    int id;
    std::vector<std::string> matches;
    if ( match( name, matches, id ) ) {
      util::Config gridconf( config );
      int CubeNx = to_int( matches[0] );
      gridconf.set( "type", type() );
      gridconf.set( "CubeNx", CubeNx );
      return create( gridconf );
    }
    return nullptr;
  }

  // Factory constructor
  const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
    int CubeNx = 0;
    config.get( "CubeNx", CubeNx );
    util::Config projconf;
    projconf.set("type", "cubedsphere_equiangular");
    projconf.set("CubeNx", CubeNx);

    // Shift projection by a longitude
    if (config.has("ShiftLon")) {
      double shiftLon = 0.0;
      config.get("ShiftLon", shiftLon);
      projconf.set("ShiftLon", shiftLon);
    }

    // Apply a Schmidt transform
    if (config.has("DoSchmidt")) {
      bool doSchmidt = false;
      config.get("DoSchmidt", doSchmidt);
      if (doSchmidt) {
        double stretchFac;
        double targetLon;
        double targetLat;
        config.get("StretchFac", stretchFac);
        config.get("TargetLon", targetLon);
        config.get("TargetLat", targetLat);
        projconf.set("DoSchmidt", doSchmidt);
        projconf.set("StretchFac", stretchFac);
        projconf.set("TargetLon", targetLon);
        projconf.set("TargetLat", targetLat);
      }
    }

    return new CubedSphereGrid::grid_t( "CS-EA-" + std::to_string( CubeNx ), CubeNx, Projection( projconf ) );
  }

  void force_link() {}

} cubedsphere_equiangular_;

// -------------------------------------------------------------------------------------------------

static class cubedsphere_equidistant : public GridBuilder {
public:
  cubedsphere_equidistant() : GridBuilder( "cubedsphere_equidistant", {"^[Cc][Ss][_-][Ee][Dd][-_]([0-9]+)$"}, {"CSED<cubedsphere>"} ) {}

  void print( std::ostream& os ) const override {
    os << std::left << std::setw( 20 ) << "CS-ED-<FaceNx>" << "Cubed sphere, equidistant";
  }

  const atlas::Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
    int id;
    std::vector<std::string> matches;
    if ( match( name, matches, id ) ) {
      util::Config gridconf( config );
      int CubeNx = to_int( matches[0] );
      gridconf.set( "type", type() );
      gridconf.set( "CubeNx", CubeNx );
      return create( gridconf );
    }
    return nullptr;
  }

  const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
    int CubeNx = 0;
    config.get( "CubeNx", CubeNx );
    util::Config projconf;
    projconf.set("type", "cubedsphere_equidistant");
    projconf.set("CubeNx", CubeNx);

    // Shift projection by a longitude
    if (config.has("ShiftLon")) {
      double shiftLon = 0.0;
      config.get("ShiftLon", shiftLon);
      projconf.set("ShiftLon", shiftLon);
    }

    // Apply a Schmidt transform
    if (config.has("DoSchmidt")) {
      bool doSchmidt = false;
      config.get("DoSchmidt", doSchmidt);
      if (doSchmidt) {
        double stretchFac;
        double targetLon;
        double targetLat;
        config.get("StretchFac", stretchFac);
        config.get("TargetLon", targetLon);
        config.get("TargetLat", targetLat);
        projconf.set("DoSchmidt", doSchmidt);
        projconf.set("StretchFac", stretchFac);
        projconf.set("TargetLon", targetLon);
        projconf.set("TargetLat", targetLat);
      }
    }

    return new CubedSphereGrid::grid_t( "CS-ED-" + std::to_string( CubeNx ), CubeNx, Projection( projconf ) );
  }

  void force_link() {}

} cubedsphere_equidistant_;

// -------------------------------------------------------------------------------------------------

static class cubedsphere_equidistant_fv3 : public GridBuilder {
public:
  cubedsphere_equidistant_fv3() : GridBuilder( "cubedsphere_equidistant_fv3", {"^[Cc][Ss][_-][Ff][Vv]3[-_]([0-9]+)$"}, {"CSFV3<cubedsphere>"} ) {}

  void print( std::ostream& os ) const override {
    os << std::left << std::setw( 20 ) << "CS-FV3-<FaceNx>" << "Cubed sphere, equidistant FV3";
  }

  const atlas::Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
    int id;
    std::vector<std::string> matches;
    if ( match( name, matches, id ) ) {
      util::Config gridconf( config );
      int CubeNx = to_int( matches[0] );
      gridconf.set( "type", type() );
      gridconf.set( "CubeNx", CubeNx );
      return create( gridconf );
    }
    return nullptr;
  }

  const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
    int CubeNx = 0;
    config.get( "CubeNx", CubeNx );
    util::Config projconf;
    projconf.set("type", "cubedsphere_equidistant_fv3");
    projconf.set("CubeNx", CubeNx);

    // Shift projection by a longitude
    double shiftLon = 0.0;
    if (config.has("ShiftLon")) {
      config.get("ShiftLon", shiftLon);
      projconf.set("ShiftLon", shiftLon);
    }

    // Apply a Schmidt transform
    bool doSchmidt = false;
    if (config.has("DoSchmidt")) {
      config.get("DoSchmidt", doSchmidt);
      if (doSchmidt) {
        double stretchFac;
        double targetLon;
        double targetLat;
        config.get("StretchFac", stretchFac);
        config.get("TargetLon", targetLon);
        config.get("TargetLat", targetLat);
        projconf.set("DoSchmidt", doSchmidt);
        projconf.set("StretchFac", stretchFac);
        projconf.set("TargetLon", targetLon);
        projconf.set("TargetLat", targetLat);
      }
    }

    return new CubedSphereGrid::grid_t( "CS-FV3-" + std::to_string( CubeNx ), CubeNx, Projection( projconf ) );
  }

  void force_link() {}

} cubedsphere_equidistant_fv3_;

// -------------------------------------------------------------------------------------------------

void force_link_CubedSphere() {
    cubedsphere_equiangular_.force_link();
    cubedsphere_equidistant_.force_link();
    cubedsphere_equidistant_fv3_.force_link();
  }

// -------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
