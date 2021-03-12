/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/interpolation/method/knn/GridBox.h"

#include <algorithm>
#include <ostream>
#include <vector>

#include "eckit/types/FloatCompare.h"
#include "eckit/types/Fraction.h"

#include "atlas/grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Earth.h"
#include "atlas/util/GaussianLatitudes.h"
#include "atlas/util/Point.h"


static constexpr double GLOBE      = 360;
static constexpr double NORTH_POLE = 90;
static constexpr double SOUTH_POLE = -90;
static constexpr double TWO        = 2;

double normalise( double lon, double minimum ) {
    while ( lon < minimum ) {
        lon += GLOBE;
    }
    while ( lon >= minimum + GLOBE ) {
        lon -= GLOBE;
    }
    return lon;
}

#define ASSERT_BOX( n, w, s, e )                                  \
    ATLAS_ASSERT( SOUTH_POLE <= s && s <= n && n <= NORTH_POLE ); \
    ATLAS_ASSERT( w <= e && e <= w + GLOBE )


namespace atlas {
namespace interpolation {
namespace method {


GridBox::GridBox( double north, double west, double south, double east ) :
    north_( north ), west_( west ), south_( south ), east_( east ) {
    ASSERT_BOX( north_, west_, south_, east_ );
}


double GridBox::area() const {
    return util::Earth::area( {west_, north_}, {east_, south_} );
}


double GridBox::diagonal() const {
    return util::Earth::distance( {west_, north_}, {east_, south_} );
}


bool GridBox::intersects( GridBox& other ) const {
    double n = std::min( north_, other.north_ );
    double s = std::max( south_, other.south_ );

    if ( !eckit::types::is_strictly_greater( n, s ) ) {
        return false;
    }

    auto intersect = []( const GridBox& a, const GridBox& b, double& w, double& e ) {
        double ref = normalise( b.west_, a.west_ );
        double w_  = std::max( a.west_, ref );
        double e_  = std::min( a.east_, normalise( b.east_, ref ) );

        if ( eckit::types::is_strictly_greater( e_, w_ ) ) {
            w = w_;
            e = e_;
            return true;
        }
        return false;
    };

    double w = std::min( west_, other.west_ );
    double e = w;

    if ( west_ <= other.west_ ? intersect( *this, other, w, e ) || intersect( other, *this, w, e )
                              : intersect( other, *this, w, e ) || intersect( *this, other, w, e ) ) {
        ATLAS_ASSERT( w <= e );
        other = {n, w, s, e};
        return true;
    }
    return false;
}


void GridBox::print( std::ostream& out ) const {
    out << "GridBox[north=" << north_ << ",west=" << west_ << ",south=" << south_ << ",east=" << east_ << "]";
}


GridBoxes::GridBoxes( const Grid& grid, bool gaussianWeightedLatitudes ) {
    StructuredGrid structured( grid );
    if ( !structured || grid.projection() ) {
        throw_NotImplemented( "GridBoxes only support structured, unprojected/unrotated grids", Here() );
    }


    // Bounding box and periodicity
    RectangularDomain domain = structured.domain();
    ATLAS_ASSERT( domain );

    auto north = domain.ymax();
    auto south = domain.ymin();
    auto west  = domain.xmin();
    auto east  = domain.xmax();
    ASSERT_BOX( north, west, south, east );

    bool periodic( ZonalBandDomain( structured.domain() ) );


    // Calculate grid-box parallels (latitude midpoints, or accumulated Gaussian quadrature weights)
    auto& y = structured.yspace();

    std::vector<double> lat;
    lat.reserve( y.size() + 1 );

    GaussianGrid gaussian;
    if ( gaussianWeightedLatitudes && ( gaussian = grid ) ) {
        auto N = [&]() {
            auto N = gaussian.N();
            ATLAS_ASSERT( N > 0 );
            return size_t( N );
        }();

        std::vector<double> latitudes( N * 2 );
        std::vector<double> weights( N * 2 );
        util::gaussian_quadrature_npole_spole( N, latitudes.data(), weights.data() );

        std::vector<double> lat_global( N * 2 + 1 );
        auto b   = lat_global.rbegin();
        auto f   = lat_global.begin();
        *( b++ ) = SOUTH_POLE;
        *( f++ ) = NORTH_POLE;

        double wacc = -1.;
        for ( size_t j = 0; j < N; ++j, ++b, ++f ) {
            wacc += TWO * weights[j];
            double deg = util::Constants::radiansToDegrees() * std::asin( wacc );
            *b         = deg;
            *f         = -( *b );
        }
        lat_global[N] = 0.;  // (equator)

        // Grids not covering the poles need clipping
        size_t j = 0;
        for ( ; j < latitudes.size(); ++j ) {
            if ( eckit::types::is_approximately_lesser_or_equal( latitudes[j], north ) &&
                 eckit::types::is_approximately_lesser_or_equal( south, latitudes[j] ) ) {
                lat.push_back( std::min( north, lat_global[j] ) );
            }
            else if ( !lat.empty() ) {
                break;
            }
        }
        ATLAS_ASSERT( !lat.empty() );
        lat.push_back( std::max( south, lat_global[j] ) );
    }
    else {
        // non-Gaussian (but structured) grids
        lat.push_back( north );
        for ( auto b = y.begin(), a = b++; b != y.end(); a = b++ ) {
            lat.push_back( ( *b + *a ) / TWO );
        }
        lat.push_back( south );
    }


    // Calculate grid-box meridians (longitude midpoints)
    auto& x = structured.xspace();
    ATLAS_ASSERT( x.nx().size() == x.dx().size() );
    ATLAS_ASSERT( x.nx().size() == x.xmin().size() );

    auto gridSize = [&]() {
        auto N = grid.size();
        ATLAS_ASSERT( N >= 0 );
        return size_t( N );
    }();

    clear();
    reserve( gridSize );
    for ( size_t j = 0; j < x.nx().size(); ++j ) {
        eckit::Fraction dx( x.dx()[j] );
        eckit::Fraction xmin( x.xmin()[j] );

        auto n = ( xmin / dx ).integralPart();
        if ( n * dx < xmin ) {
            n += 1;  // (adjust double-fraction conversions)
        }

        // On non-periodic grids, West- and East-most grid-boxes need clipping
        ATLAS_ASSERT( x.nx()[j] > 0 );

        eckit::Fraction lon0 = ( n * dx ) - ( dx / 2 );
        eckit::Fraction lon1 = lon0;
        for ( idx_t i = 0; i < x.nx()[j]; ++i ) {
            double lon0 = lon1;
            lon1 += dx;

            double n = lat[j];
            double s = lat[j + 1];
            double w = periodic ? lon0 : std::max<double>( west, lon0 );
            double e = periodic ? lon1 : std::min<double>( east, lon1 );
            emplace_back( GridBox( n, w, s, e ) );
        }

        if ( periodic ) {
            ATLAS_ASSERT( lon1 == lon0 + GLOBE );
        }
    }

    ATLAS_ASSERT( size() == gridSize );
}


GridBoxes::GridBoxes() = default;


double GridBoxes::getLongestGridBoxDiagonal() const {
    ATLAS_ASSERT( !empty() );

    double R = 0.;
    for ( auto& box : *this ) {
        R = std::max( R, box.diagonal() );
    }

    ATLAS_ASSERT( R > 0. );
    return R;
}


}  // namespace method
}  // namespace interpolation
}  // namespace atlas
