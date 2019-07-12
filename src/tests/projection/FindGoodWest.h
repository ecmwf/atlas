/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <sstream>
#include <vector>
#include "atlas/grid.h"
#include "atlas/projection/Projection.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/NormaliseLongitude.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace test {

struct WestFinder {
    // Given the corners of a grid A,B,C,D projected to lon,lat,
    //   it may happen that because of a bad normalisation the grid is
    //   inverted or crossed. This class can find a good "west" to use
    //   for normalisation. Vertical line marks a good "west" candidate.
    //   At the moment only west=0 or west=-180 are checked
    //
    //   |     A----<----D
    //   |    /         /             OK,  all angles counterclockwise
    //   |   B---->----C
    //  west
    //
    //     D---|->---A
    //    /    |   /                  NOT OK, all angles clockwise
    //   C---<-|--B
    //        west
    //
    //           |   A----<-----D
    //    / ...> |  /      >.../      NOT OK, 1 angle clockwise
    //   C---<---|-B
    //          west

    PointLonLat A, B, C, D;
    double west;
    StructuredGrid grid;
    int status;
    WestFinder( const StructuredGrid& _grid ) : grid{_grid} {
        ATLAS_ASSERT( grid, "grid must be a StructuredGrid" );
        if ( grid.y().front() > grid.y().back() ) {  // scanning mode
            A = grid.lonlat( 0, 0 );
            B = grid.lonlat( 0, grid.ny() - 1 );
            C = grid.lonlat( grid.nx( grid.ny() - 1 ) - 1, grid.ny() - 1 );
            D = grid.lonlat( grid.nx( 0 ) - 1, 0 );
        }
        else {
            A = grid.lonlat( 0, grid.ny() - 1 );
            B = grid.lonlat( 0, 0 );
            C = grid.lonlat( grid.nx( 0 ) - 1, 0 );
            D = grid.lonlat( grid.nx( grid.ny() - 1 ) - 1, grid.ny() - 1 );
        }
        if ( !find_west( 0. ) ) {
            if ( !find_west( -180. ) ) {
                if ( status == 3 ) {
                    ATLAS_THROW_EXCEPTION( "Not all grid points fall in the interval (0,360) or (-180,180):" );
                }
                else {
                    ATLAS_THROW_EXCEPTION(
                        "Could not find west that makes a valid quadrilateral from the grid corners: \n"
                        << *this );
                }
                // Perhaps we can throw more at it
            }
        }
    }
    void normalise( double _west ) {}

    bool find_west( double _west ) {
        west = _west;
        NormaliseLongitude normalised( west );
        A.lon()                    = normalised( A.lon() );
        B.lon()                    = normalised( B.lon() );
        C.lon()                    = normalised( C.lon() );
        D.lon()                    = normalised( D.lon() );
        double lon_min             = std::min( A.lon(), std::min( B.lon(), std::min( C.lon(), D.lon() ) ) );
        double lon_max             = std::max( A.lon(), std::max( B.lon(), std::max( C.lon(), D.lon() ) ) );
        double dlon                = lon_max - lon_min;
        double lon_max_with_buffer = lon_max + 0.1 * dlon;
        double lon_min_with_buffer = lon_min - 0.1 * dlon;
        double DAB                 = angle( D, A, B );
        double ABC                 = angle( A, B, C );
        double BCD                 = angle( B, C, D );
        double CDA                 = angle( C, D, A );
        if ( DAB < 0 && ABC < 0 && BCD < 0 && CDA < 0 ) {
            // angles are counterclockwise. First check OK --> Deeper check now
            if ( grid ) {
                auto in_bounds = [&]( idx_t i, idx_t j ) {
                    double lon = normalised( grid.lonlat( i, j ).lon() );
                    if ( lon > lon_max_with_buffer || lon < lon_min_with_buffer ) {
                        status = 3;
                        return false;
                    }
                    return true;
                };
                // top,bottom
                {
                    auto top_bottom = std::array<idx_t, 2>{0, grid.ny() - 1};
                    for ( idx_t j : top_bottom ) {
                        for ( idx_t i = 0; i < grid.nx( j ); ++i ) {
                            if ( !in_bounds( i, j ) )
                                return false;
                        }
                    }
                }
                // left, right
                {
                    for ( idx_t j = 0; j < grid.ny(); ++j ) {
                        auto left_right = std::array<idx_t, 2>{0, grid.nx( j ) - 1};
                        for ( idx_t i : left_right ) {
                            if ( !in_bounds( i, j ) )
                                return false;
                        }
                    }
                }
            }
            ATLAS_DEBUG( "all ok, west=" << west );
            return true;
        }
        else if ( DAB > 0 && ABC > 0 && BCD > 0 && CDA > 0 ) {
            status = 2;
            ATLAS_DEBUG( "swapped, west=" << west );
        }
        else {
            status = 1;
            ATLAS_DEBUG( "crossed, west=" << west );
        }
        return false;
    }
    static double angle( const PointLonLat& p1, const PointLonLat& p2, const PointLonLat& p3 ) {
        //     p1  p3
        //      \^/        negative = clockwise = expected
        //      p2
        PointLonLat p12 = {p1.lon() - p2.lon(), p1.lat() - p2.lat()};
        PointLonLat p32 = {p3.lon() - p2.lon(), p3.lat() - p2.lat()};
        p12 *= util::Constants::degreesToRadians();
        p32 *= util::Constants::degreesToRadians();
        double alpha =
            std::atan2( p12.lon() * p32.lat() - p12.lat() * p32.lon(), p12.lon() * p32.lon() + p12.lat() * p32.lat() ) *
            util::Constants::radiansToDegrees();
        return alpha;
    }
    void print( std::ostream& out ) const {
        out << "  A = " << A << std::endl;
        out << "  B = " << B << std::endl;
        out << "  C = " << C << std::endl;
        out << "  D = " << D << std::endl;
        out << "  DAB = " << angle( D, A, B ) << std::endl;
        out << "  ABC = " << angle( A, B, C ) << std::endl;
        out << "  BCD = " << angle( B, C, D ) << std::endl;
        out << "  CDA = " << angle( C, D, A ) << std::endl;
    }
    friend std::ostream& operator<<( std::ostream& out, const WestFinder& inst ) {
        inst.print( out );
        return out;
    }
};

//-----------------------------------------------------------------------------

double find_good_west( const Grid& grid ) {
    return WestFinder( grid ).west;
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas
