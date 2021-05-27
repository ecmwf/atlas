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
#include "atlas/projection/detail/CubedSphereProjectionBase.h"
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

CubedSphere::CubedSphere( const std::string& name, int N, Projection projection ) :
                           Grid(), N_(N), name_( name )  { // Number of tiles hardwired to 6 at the moment. Regional may need 1
  // Copy members
  projection_ = projection ? projection : Projection();

  // Domain
  domain_ = computeDomain();

  // x and y are the position in a 2D plane for the unfolded cubed-sphere grid, shown in the
  // comments in grid/CubedSphereGrid.h. In order to locate the position in this xy array the start
  // position for each face (tile) of the cube is needed. xs represents the x start position and ys
  // the y start position. Tile 3, 4 and 5 are rotated and ysr provides the start point for y after
  // these rotations.

  xs_[0] = 0*N;
  xs_[1] = 1*N;
  xs_[2] = 1*N;
  xs_[3] = 2*N;
  xs_[4] = 3*N;
  xs_[5] = 3*N;

  ys_[0] = 1*N;
  ys_[1] = 1*N;
  ys_[2] = 2*N;
  ys_[3] = 1*N+1;
  ys_[4] = 1*N+1;
  ys_[5] = 0*N+1;

  ysr_[0] = ys_[0];
  ysr_[1] = ys_[1];
  ysr_[2] = ys_[2];
  ysr_[3] = 2*N;
  ysr_[4] = 2*N;
  ysr_[5] = 1*N;

  // Number of grid points on each face of the tile.
  npts_.push_back(N*N+1); // An extra point lies on tile 1
  npts_.push_back(N*N+1); // An extra point lies on tile 2
  npts_.push_back(N*N);
  npts_.push_back(N*N);
  npts_.push_back(N*N);
  npts_.push_back(N*N);
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
  h.add(N_);

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

// Convert from xy space into resolution dependent xyt space.
void CubedSphere::xy2xyt(const double xy[], double xyt[]) const {
    // xy is in degrees while xyt is in radians
    // (alpha, beta) and tiles.

    double normalisedX = xy[XX]/90.;
    double normalisedY = (xy[YY] + 135.)/90.;

    double NDouble = static_cast<double>(N_);

    std::vector<double> yOffset{NDouble,
                                NDouble,
                                2. *  NDouble,
                                NDouble,
                                NDouble,
                                0};

    xyt[0] = (normalisedX - std::floor(normalisedX)) * static_cast<double>(N_)
          + xs_[static_cast<size_t>(xyt[2])];

    xyt[1] = (normalisedY - std::floor(normalisedY)) * static_cast<double>(N_)
          + yOffset[static_cast<size_t>(xyt[2])];

    using atlas::projection::detail::CubedSphereProjectionBase;
    xyt[2] =
        dynamic_cast<const CubedSphereProjectionBase &>(projection_).tileFromXY(xy);
}

// Convert from xyt space into continuous xy space.
void CubedSphere::xyt2xy(const double xyt[], double xy[]) const {
    // xy is in degrees
    // while xyt is in number of grid points
    // (alpha, beta) and tiles.
    std::vector<double> xOffsetDeg{0., 90., 90., 180., 270., 270.};
    std::vector<double> yOffsetDeg{-45., -45., 45., -45., -45., -135.};

    double N = static_cast<double>(N_);
    std::vector<double> xOffsetIndex{0, N, N, 2*N, 3*N,  3*N};
    std::vector<double> yOffsetIndex{N, N, 2*N, N,  N, 0};

    double normalisedX =
     (xyt[0] - xOffsetIndex[static_cast<size_t>(xyt[2])])/N;
    double normalisedY =
     (xyt[1] - yOffsetIndex[static_cast<size_t>(xyt[2])])/N;
    xy[XX] = normalisedX * 90. + xOffsetDeg[xyt[2]];
    xy[YY] = normalisedY * 90. + yOffsetDeg[xyt[2]];
}

// ------------------------------------------

namespace {
  GridFactoryBuilder<CubedSphere> __register_CubedSphere( CubedSphere::static_type() );
}

// -------------------------------------------------------------------------------------------------

// Specialization based on type of projection
// ------------------------------------------

static class cubedsphere_equiangular : public GridBuilder {
public:
  cubedsphere_equiangular() : GridBuilder( "cubedsphere_equiangular", {"^[Cc][Ss][_-][Ee][Aa][-_]([0-9]+)$"}, {"CS-EA-<N>"} ) {}

  void print( std::ostream& os ) const override {
    os << std::left << std::setw( 20 ) << "CS-EA-<N>" << "Cubed sphere, equiangular";
  }

  // Factory constructor
  const atlas::Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
    int id;
    std::vector<std::string> matches;
    if ( match( name, matches, id ) ) {
      util::Config gridconf( config );
      int N = to_int( matches[0] );
      gridconf.set( "type", type() );
      gridconf.set( "N", N );
      return create( gridconf );
    }
    return nullptr;
  }

  // Factory constructor
  const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
    int N = 0;
    if( not config.get( "N", N ) ) {
        throw_AssertionFailed("Could not find \"N\" in configuration of cubed sphere grid",Here());
    }
    util::Config projconf;
    projconf.set("type", "cubedsphere_equiangular");
    projconf.set("N", N);

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
    return new CubedSphereGrid::grid_t( "CS-EA-" + std::to_string( N ), N, Projection( projconf ) );
  }

  void force_link() {}

} cubedsphere_equiangular_;

// -------------------------------------------------------------------------------------------------

static class cubedsphere_equidistant : public GridBuilder {
public:
  cubedsphere_equidistant() : GridBuilder( "cubedsphere_equidistant", {"^[Cc][Ss][_-][Ee][Dd][-_]([0-9]+)$"}, {"CS-ED-<N>"} ) {}

  void print( std::ostream& os ) const override {
    os << std::left << std::setw( 20 ) << "CS-ED-<N>" << "Cubed sphere, equidistant";
  }

  const atlas::Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
    int id;
    std::vector<std::string> matches;
    if ( match( name, matches, id ) ) {
      util::Config gridconf( config );
      int N = to_int( matches[0] );
      gridconf.set( "type", type() );
      gridconf.set( "N", N );
      return create( gridconf );
    }
    return nullptr;
  }

  const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
    int N = 0;
    if( not config.get( "N", N ) ) {
        throw_AssertionFailed("Could not find \"N\" in configuration of cubed sphere grid",Here());
    }
    util::Config projconf;
    projconf.set("type", "cubedsphere_equidistant");
    projconf.set("N", N);

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

    return new CubedSphereGrid::grid_t( "CS-ED-" + std::to_string( N ), N, Projection( projconf ) );
  }

  void force_link() {}

} cubedsphere_equidistant_;

// -------------------------------------------------------------------------------------------------

void force_link_CubedSphere() {
    cubedsphere_equiangular_.force_link();
    cubedsphere_equidistant_.force_link();
}

Grid::Config CubedSphere::meshgenerator() const {
    return Config( "type", "cubedsphere" );
}

Grid::Config CubedSphere::partitioner() const {
    // TODO: implement better one specific for cubed sphere
    Grid::Config config;
    config.set("type","equal_regions");
    config.set("coordinates","lonlat"); // do not use the grid.xy() coordinates for partitioning
    return config;
}

// -------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
