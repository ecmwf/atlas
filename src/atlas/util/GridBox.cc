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


#include "atlas/util/GridBox.h"

#include <algorithm>
#include <ostream>

#include "eckit/types/FloatCompare.h"
#include "eckit/types/Fraction.h"

#include "atlas/grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Earth.h"
#include "atlas/util/Point.h"


static constexpr double GLOBE      = 360.;
static constexpr double NORTH_POLE = 90.;
static constexpr double SOUTH_POLE = -90.;

double normalise( double lon, double minimum ) {
    while ( lon < minimum ) {
        lon += GLOBE;
    }
    while ( lon >= minimum + GLOBE ) {
        lon -= GLOBE;
    }
    return lon;
}


namespace atlas {
namespace util {


GridBox::GridBox( double north, double west, double south, double east ) :
    north_( north ),
    west_( west ),
    south_( south ),
    east_( east ) {
    ATLAS_ASSERT( SOUTH_POLE <= south_ && south_ <= north_ && north_ <= NORTH_POLE );
    ATLAS_ASSERT( west_ <= east_ && east_ <= west_ + GLOBE );
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


GridBoxes::GridBoxes( const Grid& grid ) {
    StructuredGrid structured( grid );
    if ( !structured || grid.projection() ) {
        throw_NotImplemented( "GridBoxes only support structured, unprojected/unrotated grids", Here() );
    }


    // Calculate grid-box parallels (latitude midpoints)
    auto& y = structured.yspace();

    std::vector<double> lat;
    lat.reserve( y.size() + 1 );

    RectangularDomain domain = structured.domain();
    ATLAS_ASSERT( domain );

    auto north = domain.ymax();
    auto south = domain.ymin();
    ATLAS_ASSERT( -90. <= south && south <= north && north <= 90. );

    lat.push_back( north );
    for ( auto b = y.begin(), a = b++; b != y.end(); a = b++ ) {
        lat.push_back( ( *b + *a ) / 2. );
    }
    lat.push_back( south );

    lat.front() = std::min( north, std::max( south, lat.front() ) );  // clip to domain
    lat.back()  = std::min( north, std::max( south, lat.back() ) );   // (...)


    // Calculate grid-box meridians (longitude midpoints)
    auto& x = structured.xspace();
    ATLAS_ASSERT( x.nx().size() == x.dx().size() );
    ATLAS_ASSERT( x.nx().size() == x.xmin().size() );

    bool periodic( ZonalBandDomain( structured.domain() ) );

    clear();
    reserve( grid.size() );
    for ( size_t j = 0; j < x.nx().size(); ++j ) {
        eckit::Fraction dx( x.dx()[j] );
        eckit::Fraction xmin( x.xmin()[j] );

        auto n = ( xmin / dx ).integralPart();
        if ( n * dx < xmin ) {
            n += 1;  // (adjust double-fraction conversions)
        }

        eckit::Fraction lon1 = ( n * dx ) - ( dx / 2 );
        for ( idx_t i = 0; i < x.nx()[j]; ++i ) {
            double lon0 = lon1;
            lon1 += dx;
            emplace_back( GridBox( lat[j], lon0, lat[j + 1], lon1 ) );
        }

        if ( periodic ) {
            ATLAS_ASSERT( lon1 == xmin - ( dx / 2 ) + 360 );
        }
    }

    ATLAS_ASSERT( idx_t( size() ) == grid.size() );
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


}  // namespace util
}  // namespace atlas
