/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


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
    Grid gridA( "O32" );
    Grid gridB( "O64" );
    Grid gridC( "O32" );


    // the integral of a constant field = 1 is Earth's surface area
    struct reference_t {
        reference_t( double value, double tolerance ) : value_( value ), tolerance_( tolerance ) {
            Log::info() << "Reference: " << value_ << ", abs. tolerance: +-" << tolerance_ << std::endl;
        }

        void check( std::string label, double value ) {
            Log::info() << label << value << std::endl;
            EXPECT( eckit::types::is_approximately_equal( value, value_, tolerance_ ) );
        }

        const double value_;
        const double tolerance_;
    } reference( util::Earth::area(), 1.e2 );

    Field fieldA( "A", array::make_datatype<double>(), array::make_shape( gridA.size() ) );
    array::make_view<double, 1>( fieldA ).assign( 1. );
    reference.check( "Integral A: ", integral( gridA, fieldA ) );


    SECTION( "matrix-based" ) {
        Field fieldB( "B", array::make_datatype<double>(), array::make_shape( gridB.size() ) );
        Field fieldC( "C", array::make_datatype<double>(), array::make_shape( gridC.size() ) );

        auto config = option::type( "grid-box-average" );
        config.set( "matrix_free", false );

        Interpolation( config, gridA, gridB ).execute( fieldA, fieldB );
        reference.check( "Integral B: ", integral( gridB, fieldB ) );

        Interpolation( config, gridB, gridC ).execute( fieldB, fieldC );
        reference.check( "Integral C: ", integral( gridC, fieldC ) );
    }


    SECTION( "matrix-free" ) {
        Field fieldB( "B", array::make_datatype<double>(), array::make_shape( gridB.size() ) );
        Field fieldC( "C", array::make_datatype<double>(), array::make_shape( gridC.size() ) );

        auto config = option::type( "grid-box-average" );
        config.set( "matrix_free", true );

        Interpolation( config, gridA, gridB ).execute( fieldA, fieldB );
        reference.check( "Integral B: ", integral( gridB, fieldB ) );

        Interpolation( config, gridB, gridC ).execute( fieldB, fieldC );
        reference.check( "Integral C: ", integral( gridC, fieldC ) );
    }
}


}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
