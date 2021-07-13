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
#include "atlas/grid/Tiles.h"
#include "atlas/projection/detail/CubedSphereProjectionBase.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Point.h"
#include "atlas/util/UnitSphere.h"

static std::string extractStagger(std::string fullString) {
    // return a substring of length 1, starting one position to the left of
    // the last "-"
    return fullString.substr(fullString.rfind("-")-1, 1);
}

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

using atlas::projection::detail::CubedSphereProjectionBase;


CubedSphere::CubedSphere( int N, Projection p ) : CubedSphere( CubedSphere::static_type(), N, p ) {}

CubedSphere::CubedSphere( const std::string& name, int N, Projection projection ) :
    Grid(), N_( N ), name_( name ), stagger_( extractStagger(name) ) {
    // Number of tiles hardwired to 6 at the moment. Regional may need 1
    // Copy members
    util::Config defaultProjConfig;
    defaultProjConfig.set("type", "cubedsphere_equiangular");
    projection_ = projection ? projection : Projection(defaultProjConfig);

    // Domain
    domain_ = computeDomain();

    // x and y are the position in a 2D plane for the unfolded cubed-sphere grid, shown in the
    // comments in grid/CubedSphereGrid.h. In order to locate the position in this xy array the start
    // position for each face (tile) of the cube is needed. xs represents the x start position and ys
    // the y start position. Tile 3, 4 and 5 are rotated and ysr provides the start point for y after
    // these rotations.

    using atlas::projection::detail::CubedSphereProjectionBase;
    cs_projection_ = dynamic_cast<CubedSphereProjectionBase*>( projection_.get() );
    if( not cs_projection_ ) {
      ATLAS_THROW_EXCEPTION( "Provided projection " << projection_.type() <<
                              " is incompatible with the CubedSphere grid type" );
    }

    tiles_ = cs_projection_->getCubedSphereTiles();

    double staggerSize = (stagger_ == "C") ? 0.5 : 0.0;
    for (std::size_t i = 0; i < nTiles_; ++i) {
      // default assumes all panels start in bottom left corner or center
      xs_[i] = tiles_.xy2abOffsets()[LON][i] * N + staggerSize;
      xsr_[i] = tiles_.xy2abOffsets()[LON][i] * N + staggerSize;
      ys_[i] = tiles_.xy2abOffsets()[LAT][i] * N + staggerSize;
      ysr_[i] = tiles_.xy2abOffsets()[LAT][i] * N + staggerSize;

      // default assumes number of points on a panel is N * N
      npts_.push_back( N * N );
    }

    // default assumes jmax_ value of N-1 on all tiles
    jmax_ = std::array<idx_t,6>  {N-1,N-1,N-1,N-1,N-1,N-1};

    if (tiles_.type() == "cubedsphere_fv3") {
      // panel 3,4,5 are reversed in that they start in top left corner
      for (std::size_t i = 3; i < nTiles_; ++i) {
        if (stagger_ == "C") {
            ysr_[i] += N - 1;
        } else {
            ys_[i] += 1;
            ysr_[i] += N;
        }
      }

      // Exceptions to N * N grid points on certain tiles
      if (stagger_ == "L") {
        npts_[0] += 1;  // An extra nodal point lies on tile 1
        npts_[1] += 1;  // An extra nodal point lies on tile 2
      }

      xtile = {[this]( int i, int j, int t ) { return this->xsPlusIndex( i, t ); },
               [this]( int i, int j, int t ) { return this->xsPlusIndex( i, t ); },
               [this]( int i, int j, int t ) { return this->xsPlusIndex( i, t ); },
               [this]( int i, int j, int t ) { return this->xsPlusIndex( j, t ); },
               [this]( int i, int j, int t ) { return this->xsPlusIndex( j, t ); },
               [this]( int i, int j, int t ) { return this->xsPlusIndex( j, t ); }};

      ytile = {[this]( int i, int j, int t ) { return this->ysPlusIndex( j, t ); },
               [this]( int i, int j, int t ) { return this->ysPlusIndex( j, t ); },
               [this]( int i, int j, int t ) { return this->ysPlusIndex( j, t ); },
               [this]( int i, int j, int t ) { return this->ysrMinusIndex( i, t ); },
               [this]( int i, int j, int t ) { return this->ysrMinusIndex( i, t ); },
               [this]( int i, int j, int t ) { return this->ysrMinusIndex( i, t ); }};

      // Exceptions to jmax_ value of N-1 on certain tiles
      if (stagger_ == "L") {
          jmax_[0] = N; // Due to extra nodal point on panel 1
      }
      for ( idx_t t = 0; t < nTiles_; ++t) {
        std::size_t rowlength =  1 + jmax_[t] - jmin_[t];
        std::vector<idx_t> imaxTile(rowlength, N-1);
        std::vector<idx_t> iminTile(rowlength, 0);
        if (stagger_ == "L") {
            // extra points
            if (t == 0) imaxTile[N] = 0;
            if (t == 1) imaxTile[0] = N;
        }
        imax_.push_back(imaxTile);
        imin_.push_back(iminTile);
      }


    } else if (tiles_.type() == "cubedsphere_lfric") {

      // panel 2, 3 starts in lower right corner initially going upwards
      xs_[2] += 1;
      xs_[3] += 1;
      xsr_[2] += N-1;
      xsr_[3] += N-1;

      // panel 5 starts in upper left corner going downwards
      if (stagger_ == "L") {
        xs_[5] += 1;
        ys_[5] += 1;
      }
      ysr_[5] += N-1;

      // Exceptions to N * N grid points on certain tiles
      if (stagger_ == "L") {
        npts_[4] = (N + 1) * (N + 1); // nodal top panel includes all edges
        npts_[5] = (N - 1) * (N - 1); // nodal bottom panel excludes all edges
      }

      xtile = {[this]( int i, int j, int t ) { return this->xsPlusIndex( i, t ); },
               [this]( int i, int j, int t ) { return this->xsPlusIndex( i, t ); },
               [this]( int i, int j, int t ) { return this->xsrMinusIndex( j, t ); },
               [this]( int i, int j, int t ) { return this->xsrMinusIndex( j, t ); },
               [this]( int i, int j, int t ) { return this->xsPlusIndex( i, t ); },
               [this]( int i, int j, int t ) { return this->xsPlusIndex( j, t ); }};

      ytile = {[this]( int i, int j, int t ) { return this->ysPlusIndex( j, t ); },
               [this]( int i, int j, int t ) { return this->ysPlusIndex( j, t ); },
               [this]( int i, int j, int t ) { return this->ysPlusIndex( i, t ); },
               [this]( int i, int j, int t ) { return this->ysPlusIndex( i, t ); },
               [this]( int i, int j, int t ) { return this->ysPlusIndex( j, t ); },
               [this]( int i, int j, int t ) { return this->ysrMinusIndex( i, t ); }};

      // Exceptions to jmax_ value of N-1 on certain tiles
      if (stagger_ == "L") {
          jmax_[4] = N;
          jmax_[5] = N-2;
      }

      for ( std::size_t t = 0; t < nTiles_; ++t) {
        std::size_t rowlength =  1 + jmax_[t] - jmin_[t];
        std::vector<idx_t> imaxTile(rowlength, N-1);
        std::vector<idx_t> iminTile(rowlength, 0);
        if (stagger_ == "L") {
            if (t == 4) { std::fill_n(imaxTile.begin(),rowlength, N); }
            if (t == 5) { std::fill_n(imaxTile.begin(),rowlength, N-2); }
        }
        imax_.push_back(imaxTile);
        imin_.push_back(iminTile);
      }
    }
}

// Provide the domain for the cubed-sphere grid, which is global.
Domain CubedSphere::computeDomain() const {
    return GlobalDomain();
}

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
    h.add( "CubedSphere" );
    h.add( N_ );

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
void CubedSphere::xy2xyt( const double xy[], double xyt[] ) const {
    // xy is in degrees while xyt is in radians
    // (alpha, beta) and tiles.

    double normalisedX = xy[XX] / 90.;
    double normalisedY = ( xy[YY] + 135. ) / 90.;

    double NDouble = static_cast<double>( N_ );

    std::array<double, 6> yOffset{NDouble, NDouble, 2. * NDouble, NDouble, NDouble, 0};

    xyt[0] =
        ( normalisedX - std::floor( normalisedX ) ) * static_cast<double>( N_ ) + xs_[static_cast<size_t>( xyt[2] )];

    xyt[1] = ( normalisedY - std::floor( normalisedY ) ) * static_cast<double>( N_ ) +
             yOffset[static_cast<size_t>( xyt[2] )];
    xyt[2] = tiles_.tileFromXY(xy);

    throw std::runtime_error("error  xy2xyt");
}

// Convert from xyt space into continuous xy space.
void CubedSphere::xyt2xy( const double xyt[], double xy[] ) const {
    // xy is in degrees
    // while xyt is in number of grid points
    // (alpha, beta) and tiles.

    double N = static_cast<double>( N_ );
    std::size_t t = static_cast<std::size_t>(xyt[2]);

    double normalisedX =
     (xyt[0] - tiles_.xy2abOffsets()[XX][t] * N)/N;
    double normalisedY =
     (xyt[1] - tiles_.xy2abOffsets()[YY][t] * N)/N;
    xy[XX] = normalisedX * 90. + tiles_.ab2xyOffsets()[LON][t];
    xy[YY] = normalisedY * 90. + tiles_.ab2xyOffsets()[LAT][t];
}

// ------------------------------------------

namespace {
GridFactoryBuilder<CubedSphere> __register_CubedSphere( CubedSphere::static_type() );
}

// -------------------------------------------------------------------------------------------------


static class cubedsphere_lfric : public GridBuilder {
public:
    cubedsphere_lfric() :
        GridBuilder( "cubedsphere_lfric", {"^[Cc][Ss][_-][Ll][Ff][Rr][-_][CL][-_]([0-9]+)$"}, {"CS-LFR-<S>-<N>"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 20 ) << "CS-LFR-<S>-<N>"
           << "Cubed sphere for equiangular";
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
            gridconf.set( "stagger", extractStagger(name));
            return create( gridconf );
        }
        return nullptr;
    }

    // Factory constructor
    const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
        int N = 0;
        if ( not config.get( "N", N ) ) {
            throw_AssertionFailed( "Could not find \"N\" in configuration of cubed sphere grid", Here() );
        }
        std::string stagger;
        if ( not config.get( "stagger", stagger ) ) {
            stagger = "L"; // Default to nodal
        }

        util::Config projconf;
        projconf.set( "type", "cubedsphere_equiangular" );
        projconf.set( "tile.type", "cubedsphere_lfric" );

        // Shift projection by a longitude
        if ( config.has( "ShiftLon" ) ) {
            double shiftLon = 0.0;
            config.get( "ShiftLon", shiftLon );
            projconf.set( "ShiftLon", shiftLon );
        }

        // Apply a Schmidt transform
        if ( config.has( "DoSchmidt" ) ) {
            bool doSchmidt = false;
            config.get( "DoSchmidt", doSchmidt );
            if ( doSchmidt ) {
                double stretchFac;
                double targetLon;
                double targetLat;
                config.get( "StretchFac", stretchFac );
                config.get( "TargetLon", targetLon );
                config.get( "TargetLat", targetLat );
                projconf.set( "DoSchmidt", doSchmidt );
                projconf.set( "StretchFac", stretchFac );
                projconf.set( "TargetLon", targetLon );
                projconf.set( "TargetLat", targetLat );
            }
        }

        return new CubedSphereGrid::grid_t( "CS-LFR-" + stagger + "-" + std::to_string( N ),
            N, Projection( projconf ) );
    }

    void force_link() {}

} cubedsphere_lfric_;





// Specialization based on type of projection
// ------------------------------------------
static class cubedsphere_equiangular : public GridBuilder {
public:
    cubedsphere_equiangular() :
        GridBuilder( "cubedsphere_equiangular", {"^[Cc][Ss][_-][Ee][Aa][-_][CL][-_]([0-9]+)$"}, {"CS-EA-<S>-<N>"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 20 ) << "CS-EA-<S>-<N>"
           << "Cubed sphere for equiangular";
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
            gridconf.set( "stagger", extractStagger(name));
            return create( gridconf );
        }
        return nullptr;
    }

    // Factory constructor
    const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
        int N = 0;
        if ( not config.get( "N", N ) ) {
            throw_AssertionFailed( "Could not find \"N\" in configuration of cubed sphere grid", Here() );
        }
        std::string stagger;
        if ( not config.get( "stagger", stagger ) ) {
            stagger = "L"; // Default to nodal
        }
        util::Config projconf;
        projconf.set( "type", "cubedsphere_equiangular" );
        projconf.set( "tile.type", "cubedsphere_fv3" );

        // Shift projection by a longitude
        if ( config.has( "ShiftLon" ) ) {
            double shiftLon = 0.0;
            config.get( "ShiftLon", shiftLon );
            projconf.set( "ShiftLon", shiftLon );
        }

        // Apply a Schmidt transform
        if ( config.has( "DoSchmidt" ) ) {
            bool doSchmidt = false;
            config.get( "DoSchmidt", doSchmidt );
            if ( doSchmidt ) {
                double stretchFac;
                double targetLon;
                double targetLat;
                config.get( "StretchFac", stretchFac );
                config.get( "TargetLon", targetLon );
                config.get( "TargetLat", targetLat );
                projconf.set( "DoSchmidt", doSchmidt );
                projconf.set( "StretchFac", stretchFac );
                projconf.set( "TargetLon", targetLon );
                projconf.set( "TargetLat", targetLat );
            }
        }
        return new CubedSphereGrid::grid_t( "CS-EA-" + stagger + "-" + std::to_string( N ),
            N, Projection( projconf ) );
    }

    void force_link() {}

} cubedsphere_equiangular_;

// -------------------------------------------------------------------------------------------------

static class cubedsphere_equidistant : public GridBuilder {
public:
    cubedsphere_equidistant() :
        GridBuilder( "cubedsphere_equidistant", {"^[Cc][Ss][_-][Ee][Dd][-_][CL][-_]([0-9]+)$"}, {"CS-ED-<S>-<N>"} ) {}

    void print( std::ostream& os ) const override {
        os << std::left << std::setw( 20 ) << "CS-ED-<S>-<N>"
           << "Cubed sphere, equidistant";
    }

    const atlas::Grid::Implementation* create( const std::string& name, const Grid::Config& config ) const override {
        int id;
        std::vector<std::string> matches;
        if ( match( name, matches, id ) ) {
            util::Config gridconf( config );
            int N = to_int( matches[0] );
            gridconf.set( "type", type() );
            gridconf.set( "N", N );
            gridconf.set( "stagger", extractStagger(name));
            return create( gridconf );
        }
        return nullptr;
    }

    const atlas::Grid::Implementation* create( const Grid::Config& config ) const override {
        int N = 0;
        if ( not config.get( "N", N ) ) {
            throw_AssertionFailed( "Could not find \"N\" in configuration of cubed sphere grid", Here() );
        }
        std::string stagger;
        if ( not config.get( "stagger", stagger ) ) {
            stagger = "L"; // Default to nodal
        }
        util::Config projconf;
        projconf.set( "type", "cubedsphere_equidistant" );
        projconf.set( "tile.type", "cubedsphere_fv3" );

        // Shift projection by a longitude
        if ( config.has( "ShiftLon" ) ) {
            double shiftLon = 0.0;
            config.get( "ShiftLon", shiftLon );
            projconf.set( "ShiftLon", shiftLon );
        }

        // Apply a Schmidt transform
        if ( config.has( "DoSchmidt" ) ) {
            bool doSchmidt = false;
            config.get( "DoSchmidt", doSchmidt );
            if ( doSchmidt ) {
                double stretchFac;
                double targetLon;
                double targetLat;
                config.get( "StretchFac", stretchFac );
                config.get( "TargetLon", targetLon );
                config.get( "TargetLat", targetLat );
                projconf.set( "DoSchmidt", doSchmidt );
                projconf.set( "StretchFac", stretchFac );
                projconf.set( "TargetLon", targetLon );
                projconf.set( "TargetLat", targetLat );
            }
        }

        return new CubedSphereGrid::grid_t( "CS-ED-" + stagger + "-" + std::to_string( N ),
            N, Projection( projconf ) );
    }

    void force_link() {}

} cubedsphere_equidistant_;

// -------------------------------------------------------------------------------------------------

void force_link_CubedSphere() {
    cubedsphere_lfric_.force_link();
    cubedsphere_equiangular_.force_link();
    cubedsphere_equidistant_.force_link();
}

Grid::Config CubedSphere::meshgenerator() const {
    return Config( "type", "cubedsphere" );
}

Grid::Config CubedSphere::partitioner() const {
    // TODO: implement better one specific for cubed sphere
    Grid::Config config;
    config.set( "type", "equal_regions" );
    config.set( "coordinates", "lonlat" );  // do not use the grid.xy() coordinates for partitioning
    return config;
}

// -------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
