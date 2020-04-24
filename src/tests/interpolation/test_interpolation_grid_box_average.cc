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


void reset_fields( Field& f1, Field& f2, Field& f3 ) {
    array::make_view<double, 1>( f1 ).assign( 0. );
    array::make_view<double, 1>( f2 ).assign( 0. );
    array::make_view<double, 1>( f3 ).assign( 0. );
}


void reset_fieldsets( FieldSet& f1, FieldSet& f2, FieldSet& f3 ) {
    ATLAS_ASSERT( f1.size() == f2.size() );
    ATLAS_ASSERT( f1.size() == f3.size() );

    double value = 1.;  // or your money back :-)
    for ( idx_t i = 0; i < f1.size(); ++i, value += 1. ) {
        array::make_view<double, 1>( f1[i] ).assign( value );
        array::make_view<double, 1>( f2[i] ).assign( 0. );
        array::make_view<double, 1>( f3[i] ).assign( 0. );
    }
}


FieldSet create_fieldset( std::string name, idx_t size, size_t number ) {
    ATLAS_ASSERT( 1 <= number );

    FieldSet set;
    for ( size_t i = 1; i <= number; ++i ) {
        set.add( Field( name + std::to_string( i ), array::make_datatype<double>(), array::make_shape( size ) ) );
    }
    return set;
}


CASE( "test_interpolation_grid_box_average" ) {
    Log::info().precision( 16 );

    Grid gridA( "O32" );
    Grid gridB( "O64" );
    Grid gridC( "O32" );


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
    };


    // the integral of a constant field = 1 is Earth's surface area (and tolerance has to be large)
    reference_t surface1( util::Earth::area(), 1.e3 );
    reference_t surface2( util::Earth::area() * 2., 1.e3 );
    reference_t surface3( util::Earth::area() * 3., 1.e3 );

    reference_t countA( double( gridA.size() ), 1.e-8 );
    reference_t countB( double( gridB.size() ), 1.e-8 );


    // setup fields
    Field fieldA( "A", array::make_datatype<double>(), array::make_shape( gridA.size() ) );
    Field fieldB( "B", array::make_datatype<double>(), array::make_shape( gridB.size() ) );
    Field fieldC( "C", array::make_datatype<double>(), array::make_shape( gridC.size() ) );

    FieldSet fieldsA( create_fieldset( "A", gridA.size(), 3 ) );
    FieldSet fieldsB( create_fieldset( "B", gridB.size(), 3 ) );
    FieldSet fieldsC( create_fieldset( "C", gridC.size(), 3 ) );


    SECTION( "Earth's surface area, Field interpolation, matrix-based" ) {
        reset_fields( fieldA, fieldB, fieldC );
        array::make_view<double, 1>( fieldA ).assign( 1. );

        surface1.check( "Integral A: ", integral( gridA, fieldA ) );  // checked once

        auto config = option::type( "grid-box-average" );
        config.set( "matrix_free", false );

        Interpolation( config, gridA, gridB ).execute( fieldA, fieldB );
        surface1.check( "Integral B: ", integral( gridB, fieldB ) );

        Interpolation( config, gridB, gridC ).execute( fieldB, fieldC );
        surface1.check( "Integral C: ", integral( gridC, fieldC ) );
    }


    SECTION( "Earth's surface area, FieldSet interpolation, matrix-based" ) {
        reset_fieldsets( fieldsA, fieldsB, fieldsC );

        auto config = option::type( "grid-box-average" );
        config.set( "matrix_free", false );

        Interpolation( config, gridA, gridB ).execute( fieldsA, fieldsB );
        surface1.check( "Integral B1: ", integral( gridB, fieldsB[0] ) );
        surface2.check( "Integral B2: ", integral( gridB, fieldsB[1] ) );
        surface3.check( "Integral B3: ", integral( gridB, fieldsB[2] ) );

        Interpolation( config, gridB, gridC ).execute( fieldsB, fieldsC );
        surface1.check( "Integral C1: ", integral( gridC, fieldsC[0] ) );
        surface2.check( "Integral C2: ", integral( gridC, fieldsC[1] ) );
        surface3.check( "Integral C3: ", integral( gridC, fieldsC[2] ) );
    }


    SECTION( "Earth's surface area, Field interpolation, matrix-free" ) {
        reset_fields( fieldA, fieldB, fieldC );
        array::make_view<double, 1>( fieldA ).assign( 1. );

        auto config = option::type( "grid-box-average" );
        config.set( "matrix_free", true );

        Interpolation( config, gridA, gridB ).execute( fieldA, fieldB );
        surface1.check( "Integral B: ", integral( gridB, fieldB ) );

        Interpolation( config, gridB, gridC ).execute( fieldB, fieldC );
        surface1.check( "Integral C: ", integral( gridC, fieldC ) );
    }


    SECTION( "Earth's surface area, FieldSet interpolation, matrix-free" ) {
        reset_fieldsets( fieldsA, fieldsB, fieldsC );

        auto config = option::type( "grid-box-average" );
        config.set( "matrix_free", true );

        Interpolation( config, gridA, gridB ).execute( fieldsA, fieldsB );
        surface1.check( "Integral B1: ", integral( gridB, fieldsB[0] ) );
        surface2.check( "Integral B2: ", integral( gridB, fieldsB[1] ) );
        surface3.check( "Integral B3: ", integral( gridB, fieldsB[2] ) );

        Interpolation( config, gridB, gridC ).execute( fieldsB, fieldsC );
        surface1.check( "Integral C1: ", integral( gridC, fieldsC[0] ) );
        surface2.check( "Integral C2: ", integral( gridC, fieldsC[1] ) );
        surface3.check( "Integral C3: ", integral( gridC, fieldsC[2] ) );
    }


    SECTION( "Conserve count as integral value, low >> high >> low resolution, matrix-free" ) {
        reset_fields( fieldA, fieldB, fieldC );

        size_t i    = 0;
        auto values = array::make_view<double, 1>( fieldA );
        for ( auto& box : util::GridBoxes( gridA ) ) {
            ATLAS_ASSERT( box.area() > 0. );
            values( i++ ) = 1. / box.area();
        }

        auto config = option::type( "grid-box-average" );
        config.set( "matrix_free", true );

        countA.check( "Count A: ", integral( gridA, fieldA ) );

        Interpolation( config, gridA, gridB ).execute( fieldA, fieldB );
        countA.check( "Count B: ", integral( gridB, fieldB ) );

        Interpolation( config, gridB, gridA ).execute( fieldB, fieldC );
        countA.check( "Count A: ", integral( gridA, fieldC ) );
    }


    SECTION( "Conserve count as integral value, high >> low >> high resolution, matrix-free" ) {
        reset_fields( fieldA, fieldB, fieldC );

        Field fieldD( "D", array::make_datatype<double>(), array::make_shape( gridB.size() ) );
        array::make_view<double, 1>( fieldD ).assign( 0. );

        size_t i    = 0;
        auto values = array::make_view<double, 1>( fieldB );
        for ( auto& box : util::GridBoxes( gridB ) ) {
            ATLAS_ASSERT( box.area() > 0. );
            values( i++ ) = 1. / box.area();
        }

        auto config = option::type( "grid-box-average" );
        config.set( "matrix_free", true );

        countB.check( "Count B: ", integral( gridB, fieldB ) );

        Interpolation( config, gridB, gridA ).execute( fieldB, fieldA );
        countB.check( "Count A: ", integral( gridA, fieldA ) );

        Interpolation( config, gridA, gridB ).execute( fieldA, fieldD );
        countB.check( "Count B: ", integral( gridB, fieldD ) );
    }
}


}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
