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


#include <limits>

#include "atlas/array.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/interpolation/NonLinear.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/runtime/Exception.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {


const double missingValue    = 42.;
const double missingValueEps = 1e-9;
const double nan             = std::numeric_limits<double>::quiet_NaN();


CASE( "test_interpolation_non_linear_missing_value" ) {
    using interpolation::MissingValue;


    SECTION( "nan" ) {
        util::Config config;
        auto mv_nan = MissingValue( "nan", config );
        EXPECT( bool( mv_nan ) );

        EXPECT( mv_nan( nan ) );
        EXPECT( mv_nan( missingValue ) == false );

        config.set( "type", "nan" );
        mv_nan = MissingValue( config );
        EXPECT( bool( mv_nan ) );

        EXPECT( mv_nan( nan ) );
        EXPECT( mv_nan( missingValue ) == false );
    }


    SECTION( "equals" ) {
        util::Config config;
        config.set( "missing_value", missingValue );

        auto mv = MissingValue( "equals", config );
        EXPECT( bool( mv ) );

        EXPECT( mv( missingValue - 1 ) == false );
        EXPECT( mv( missingValue + 1 ) == false );
        EXPECT( mv( missingValue ) );

        config.set( "type", "equals" );
        mv = MissingValue( config );
        EXPECT( bool( mv ) );

        EXPECT( mv( missingValue - 1 ) == false );
        EXPECT( mv( missingValue + 1 ) == false );
        EXPECT( mv( missingValue ) );
    }


    SECTION( "approximately-equals" ) {
        util::Config config;
        config.set( "missing_value", missingValue );
        config.set( "missing_value_epsilon", missingValueEps );

        auto mv = MissingValue( "approximately-equals", config );
        EXPECT( bool( mv ) );

        EXPECT( mv( missingValue - missingValueEps * 2 ) == false );
        EXPECT( mv( missingValue - missingValueEps / 2 ) );
        EXPECT( mv( missingValue ) );
        EXPECT( mv( missingValue + missingValueEps / 2 ) );
        EXPECT( mv( missingValue + missingValueEps * 2 ) == false );

        config.set( "type", "approximately-equals" );
        mv = MissingValue( config );
        EXPECT( bool( mv ) );

        EXPECT( mv( missingValue - missingValueEps * 2 ) == false );
        EXPECT( mv( missingValue - missingValueEps / 2 ) );
        EXPECT( mv( missingValue ) );
        EXPECT( mv( missingValue + missingValueEps / 2 ) );
        EXPECT( mv( missingValue + missingValueEps * 2 ) == false );
    }
}


CASE( "test_interpolation_non_linear_matrix" ) {
    /*
       Set input field full of 1's, with 9 nodes
         1 ... 1 ... 1
         :     :     :
         1-----m ... 1  m: missing value
         |i   i|     :  i: interpolation on two points, this quadrilateral only
         1-----1 ... 1
     */
    RectangularDomain domain( {0, 2}, {0, 2}, "degrees" );
    Grid gridA( "L90", domain );

    const idx_t nbNodes = 9;
    ATLAS_ASSERT( gridA.size() == nbNodes );

    Mesh meshA = MeshGenerator( "structured" ).generate( gridA );

    functionspace::NodeColumns fsA( meshA );
    Field fieldA = fsA.createField<double>( option::name( "A" ) );

    auto viewA = array::make_view<double, 1>( fieldA );
    for ( idx_t j = 0; j < fsA.nodes().size(); ++j ) {
        viewA( j ) = 1;
    }


    // Set output field (2 points)
    functionspace::PointCloud fsB( {PointLonLat{0.1, 0.1}, PointLonLat{0.9, 0.9}} );
    Field fieldB( "B", array::make_datatype<double>(), array::make_shape( fsB.size() ) );

    auto viewB = array::make_view<double, 1>( fieldB );
    ATLAS_ASSERT( viewB.size() == 2 );


    util::Config cfg( "type", "finite-element" );
    cfg.set( "missing_value", missingValue );
    cfg.set( "missing_value_epsilon", missingValueEps );


    // NOTE: "equals" is not tested due to internal conversions
    SECTION( "missing-if-all-missing" ) {
        for ( std::string cmp : {"approximately-equals", "nan"} ) {
            viewA( 4 ) = cmp == "nan" ? nan : missingValue;

            cfg.set( "non_linear", "missing-if-all-missing" );
            cfg.set( "missing_value_compare", cmp );

            Interpolation interpolation( cfg, fsA, fsB );
            interpolation.execute( fieldA, fieldB );

            interpolation::MissingValue miss( cmp, cfg );
            EXPECT( miss( viewB( 0 ) ) == false );
            EXPECT( miss( viewB( 1 ) ) == false );
        }
    }


    SECTION( "missing-if-any-missing" ) {
        for ( std::string cmp : {"approximately-equals", "nan"} ) {
            viewA( 4 ) = cmp == "nan" ? nan : missingValue;

            cfg.set( "non_linear", "missing-if-any-missing" );
            cfg.set( "missing_value_compare", cmp );

            Interpolation interpolation( cfg, fsA, fsB );
            interpolation.execute( fieldA, fieldB );

            interpolation::MissingValue miss( cmp, cfg );
            EXPECT( miss( viewB( 0 ) ) );
            EXPECT( miss( viewB( 1 ) ) );
        }
    }


    SECTION( "missing-if-heaviest-missing" ) {
        for ( std::string cmp : {"approximately-equals", "nan"} ) {
            viewA( 4 ) = cmp == "nan" ? nan : missingValue;

            cfg.set( "non_linear", "missing-if-heaviest-missing" );
            cfg.set( "missing_value_compare", cmp );

            Interpolation interpolation( cfg, fsA, fsB );
            interpolation.execute( fieldA, fieldB );

            interpolation::MissingValue miss( cmp, cfg );
            EXPECT( miss( viewB( 0 ) ) == false );
            EXPECT( miss( viewB( 1 ) ) );
        }
    }
}


}  // namespace test
}  // namespace atlas


int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
