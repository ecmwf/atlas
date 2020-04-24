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

#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"
#include "atlas/util/GridBox.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {


double integral( const Grid& grid, const Field& field ) {
    auto values = array::make_view<double, 1>( field );

    auto boxes = util::GridBoxes( grid );
    ATLAS_ASSERT( boxes.size() == size_t( field.shape( 0 ) ) );

    double i = 0.;
    for ( size_t c = 0; c < boxes.size(); c++ ) {
        i += boxes[c].area() * values( c );
    }
    return i;
}


CASE( "test_interpolation_grid_box_average" ) {
    using eckit::types::is_approximately_equal;

    Grid gridA( "O32" );
    Grid gridB( "O64" );
    Grid gridC( "O32" );


    // the integral of field = 1 is Earth's surface area
    auto ref = util::Earth::area();
    Log::info() << "Earth's surface = " << ref << " m2 (reference)" << std::endl;
    static double interpolation_tolerance = 1e3;  // (Earth's surface area is large)

    Field fieldA( "A", array::make_datatype<double>(), array::make_shape( gridA.size() ) );
    array::make_view<double, 1>( fieldA ).assign( 1. );


    SECTION( "matrix-free" ) {
        Field fieldB( "B", array::make_datatype<double>(), array::make_shape( gridB.size() ) );
        Field fieldC( "C", array::make_datatype<double>(), array::make_shape( gridC.size() ) );

        auto config = option::type( "grid-box-average" );
        config.set( "matrix_free", true );

        Interpolation( config, gridA, gridB ).execute( fieldA, fieldB );
        Interpolation( config, gridB, gridC ).execute( fieldB, fieldC );

        double integrals[]{integral( gridA, fieldA ), integral( gridB, fieldB ), integral( gridC, fieldC )};
        Log::info() << "Integral A: " << integrals[0] << "\n"
                    << "Integral B: " << integrals[1] << "\n"
                    << "Integral C: " << integrals[2] << std::endl;

        EXPECT( is_approximately_equal( integrals[0], ref, interpolation_tolerance ) );
        EXPECT( is_approximately_equal( integrals[1], ref, interpolation_tolerance ) );
        EXPECT( is_approximately_equal( integrals[2], ref, interpolation_tolerance ) );
    }

    SECTION( "matrix-based" ) {
        Field fieldB( "B", array::make_datatype<double>(), array::make_shape( gridB.size() ) );
        Field fieldC( "C", array::make_datatype<double>(), array::make_shape( gridC.size() ) );

        auto config = option::type( "grid-box-average" );
        config.set( "matrix_free", false );

        Interpolation( config, gridA, gridB ).execute( fieldA, fieldB );
        Interpolation( config, gridB, gridC ).execute( fieldB, fieldC );

        double integrals[]{integral( gridA, fieldA ), integral( gridB, fieldB ), integral( gridC, fieldC )};
        Log::info() << "Integral A: " << integrals[0] << "\n"
                    << "Integral B: " << integrals[1] << "\n"
                    << "Integral C: " << integrals[2] << std::endl;

        EXPECT( is_approximately_equal( integrals[0], ref, interpolation_tolerance ) );
        EXPECT( is_approximately_equal( integrals[1], ref, interpolation_tolerance ) );
        EXPECT( is_approximately_equal( integrals[2], ref, interpolation_tolerance ) );
    }
}


}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
