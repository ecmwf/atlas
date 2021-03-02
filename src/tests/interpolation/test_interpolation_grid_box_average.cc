/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include <cmath>

#include "eckit/log/Bytes.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/interpolation/method/knn/GridBox.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {


using interpolation::method::GridBoxes;


double integral( const Grid& grid, const Field& field ) {
    auto values = array::make_view<double, 1>( field );

    auto boxes = GridBoxes( grid );
    ATLAS_ASSERT( boxes.size() == size_t( field.shape( 0 ) ) );

    double i = 0.;
    for ( size_t c = 0; c < boxes.size(); c++ ) {
        i += boxes[c].area() * values( c );
    }
    return i;
}


Field create_field( std::string name, idx_t size, double init = double() ) {
    auto f = Field( name, array::make_datatype<double>(), array::make_shape( size ) );
    array::make_view<double, 1>( f ).assign( init );
    return f;
}


FieldSet create_fieldset( std::string name, idx_t size, size_t number ) {
    ATLAS_ASSERT( 1 <= number );

    FieldSet set;
    for ( size_t i = 1; i <= number; ++i ) {
        auto f =
            set.add( Field( name + std::to_string( i ), array::make_datatype<double>(), array::make_shape( size ) ) );
        array::make_view<double, 1>( f ).assign( 0. );
    }
    return set;
}


CASE( "test_interpolation_grid_box_average" ) {
    Log::info().precision( 16 );

    Grid gridA( "O32" );
    Grid gridB( "O64" );
    Grid gridC( "O32" );


    struct reference_t {
        reference_t( double value, double tolerance ) : value_( value ), tolerance_( tolerance ) {}

        void check( std::string label, double value ) {
            double relerr =
                eckit::types::is_approximately_equal( value_, 0. ) ? 0. : std::abs( ( value - value_ ) / value_ );
            Log::info() << label << value << ", compared to " << value_ << " +- " << tolerance_
                        << ", with relative error = " << relerr << std::endl;
            EXPECT( eckit::types::is_approximately_equal( value, value_, tolerance_ ) );
        }

        const double value_;
        const double tolerance_;
    };


    // the integral of a constant field = 1 is Earth's surface area (and tolerance has to be large)
    // (also field = 2 and 3 just to test a bit further, too)
    reference_t surface1( util::Earth::area(), 1.e3 );
    reference_t surface2( util::Earth::area() * 2., 1.e3 );
    reference_t surface3( util::Earth::area() * 3., 1.e3 );

    reference_t countA( double( gridA.size() ), 1.e-8 );
    reference_t countB( double( gridB.size() ), 1.e-8 );

    auto configMB = option::type( "grid-box-average" ).set( "matrix_free", false );
    auto configMF = option::type( "grid-box-average" ).set( "matrix_free", true );


    SECTION( "Earth's surface area using Field interpolation" ) {
        Field fieldA( create_field( "A", gridA.size(), 1. ) );
        Field fieldB( create_field( "B", gridB.size() ) );
        Field fieldC( create_field( "C", gridC.size() ) );

        surface1.check( "Integral A: ", integral( gridA, fieldA ) );

        Interpolation( configMB, gridA, gridB ).execute( fieldA, fieldB );
        Interpolation( configMB, gridB, gridC ).execute( fieldB, fieldC );

        surface1.check( "Integral B (MB): ", integral( gridB, fieldB ) );
        surface1.check( "Integral C (MB): ", integral( gridC, fieldC ) );

        fieldB = create_field( "B", gridB.size() );
        fieldC = create_field( "C", gridC.size() );

        Interpolation( configMF, gridA, gridB ).execute( fieldA, fieldB );
        Interpolation( configMF, gridB, gridC ).execute( fieldB, fieldC );

        surface1.check( "Integral B (MF): ", integral( gridB, fieldB ) );
        surface1.check( "Integral C (MF): ", integral( gridC, fieldC ) );
    }


    SECTION( "Earth's surface area using FieldSet interpolation" ) {
        FieldSet fieldsA( create_fieldset( "A", gridA.size(), 3 ) );
        array::make_view<double, 1>( fieldsA[0] ).assign( 1. );
        array::make_view<double, 1>( fieldsA[1] ).assign( 2. );
        array::make_view<double, 1>( fieldsA[2] ).assign( 3. );

        FieldSet fieldsB( create_fieldset( "B", gridB.size(), 3 ) );
        FieldSet fieldsC( create_fieldset( "C", gridC.size(), 3 ) );

        surface1.check( "Integral A1: ", integral( gridA, fieldsA[0] ) );
        surface2.check( "Integral A2: ", integral( gridA, fieldsA[1] ) );
        surface3.check( "Integral A3: ", integral( gridA, fieldsA[2] ) );

        Interpolation( configMB, gridA, gridB ).execute( fieldsA, fieldsB );
        Interpolation( configMB, gridB, gridC ).execute( fieldsB, fieldsC );

        surface1.check( "Integral B1 (MB): ", integral( gridB, fieldsB[0] ) );
        surface2.check( "Integral B2 (MB): ", integral( gridB, fieldsB[1] ) );
        surface3.check( "Integral B3 (MB): ", integral( gridB, fieldsB[2] ) );

        surface1.check( "Integral C1 (MB): ", integral( gridC, fieldsC[0] ) );
        surface2.check( "Integral C2 (MB): ", integral( gridC, fieldsC[1] ) );
        surface3.check( "Integral C3 (MB): ", integral( gridC, fieldsC[2] ) );

        fieldsB = create_fieldset( "B", gridB.size(), 3 );
        fieldsC = create_fieldset( "C", gridC.size(), 3 );

        Interpolation( configMF, gridA, gridB ).execute( fieldsA, fieldsB );
        Interpolation( configMF, gridB, gridC ).execute( fieldsB, fieldsC );

        surface1.check( "Integral B1 (MF): ", integral( gridB, fieldsB[0] ) );
        surface2.check( "Integral B2 (MF): ", integral( gridB, fieldsB[1] ) );
        surface3.check( "Integral B3 (MF): ", integral( gridB, fieldsB[2] ) );

        surface1.check( "Integral C1 (MF): ", integral( gridC, fieldsC[0] ) );
        surface2.check( "Integral C2 (MF): ", integral( gridC, fieldsC[1] ) );
        surface3.check( "Integral C3 (MF): ", integral( gridC, fieldsC[2] ) );
    }


    SECTION( "Count as integral value, low >> high >> low resolution" ) {
        Field fieldA( create_field( "A", gridA.size(), 1. ) );
        Field fieldB( create_field( "B", gridB.size() ) );
        Field fieldC( create_field( "C", gridA.size() ) );

        size_t i    = 0;
        auto values = array::make_view<double, 1>( fieldA );
        for ( auto& box : GridBoxes( gridA ) ) {
            ATLAS_ASSERT( box.area() > 0. );
            values( i++ ) = 1. / box.area();
        }

        countA.check( "Count A: ", integral( gridA, fieldA ) );

        Interpolation( configMB, gridA, gridB ).execute( fieldA, fieldB );
        Interpolation( configMB, gridB, gridA ).execute( fieldB, fieldC );

        countA.check( "Count B (MB): ", integral( gridB, fieldB ) );
        countA.check( "Count A (MB): ", integral( gridA, fieldC ) );
    }


    SECTION( "Count as integral value, high >> low >> high resolution" ) {
        Field fieldA( create_field( "A", gridB.size(), 1. ) );
        Field fieldB( create_field( "B", gridA.size() ) );
        Field fieldC( create_field( "C", gridB.size() ) );

        size_t i    = 0;
        auto values = array::make_view<double, 1>( fieldA );
        for ( auto& box : GridBoxes( gridB ) ) {
            ATLAS_ASSERT( box.area() > 0. );
            values( i++ ) = 1. / box.area();
        }

        countB.check( "Count B: ", integral( gridB, fieldA ) );

        Interpolation( configMB, gridB, gridA ).execute( fieldA, fieldB );
        Interpolation( configMB, gridA, gridB ).execute( fieldB, fieldC );

        countB.check( "Count A (MB): ", integral( gridA, fieldB ) );
        countB.check( "Count B (MB): ", integral( gridB, fieldC ) );

        fieldB = create_field( "B", gridA.size() );
        fieldC = create_field( "C", gridB.size() );

        Interpolation( configMF, gridB, gridA ).execute( fieldA, fieldB );
        Interpolation( configMF, gridA, gridB ).execute( fieldB, fieldC );

        countB.check( "Count A (MF): ", integral( gridA, fieldB ) );
        countB.check( "Count B (MF): ", integral( gridB, fieldC ) );
    }
}

CASE( "test_interpolation_grid_box_average matrix with cache" ) {
    Grid gridA( "O32" );
    Grid gridB( "O64" );

    auto config = option::type( "grid-box-average" ).set( "matrix_free", false );

    Field fieldA( create_field( "A", gridA.size(), 1. ) );
    Field fieldB( create_field( "B", gridB.size() ) );

    size_t i    = 0;
    auto values = array::make_view<double, 1>( fieldA );
    for ( auto& box : GridBoxes( gridA ) ) {
        ATLAS_ASSERT( box.area() > 0. );
        values( i++ ) = 1. / box.area();
    }

    interpolation::Cache cache;
    ATLAS_TRACE_SCOPE( "Create cache" ) { cache = interpolation::Cache( Interpolation( config, gridA, gridB ) ); }
    size_t size;
    size_t footprint;
    ATLAS_TRACE_SCOPE( "iterate tree" ) {
        size      = interpolation::IndexKDTreeCache( cache ).tree().size();
        footprint = interpolation::IndexKDTreeCache( cache ).tree().footprint();
    }
    ATLAS_DEBUG_VAR( size );
    Log::info() << "kdtree footprint = " << eckit::Bytes( footprint ) << std::endl;

    Log::info() << "all cache footprint = " << eckit::Bytes( cache.footprint() ) << std::endl;

    ATLAS_TRACE_SCOPE( "Interpolate with cache" ) {
        Interpolation( config, gridA, gridB, cache ).execute( fieldA, fieldB );
    }
}

}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
