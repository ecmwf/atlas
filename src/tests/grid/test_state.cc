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
#include <string>

#include "eckit/eckit_version.h"
#include "eckit/log/JSON.h"

#include "eckit/parser/JSONParser.h"

#include "atlas/array/ArrayView.h"
#include "atlas/array/DataType.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/State.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::field;
using namespace atlas::field;

namespace atlas {
namespace test {

// -------------------------------------------------------------------
// Example State generator

// ---  Declaration (in .h file)
class MyStateGenerator : public StateGenerator {
public:
    MyStateGenerator( const eckit::Parametrisation& p = util::Config() ) : StateGenerator( p ) {}
    ~MyStateGenerator() override = default;
    void generate( State& state, const eckit::Parametrisation& p = util::Config() ) const override;
};

// ---  Implementation (in .cc file)
void MyStateGenerator::generate( State& state, const eckit::Parametrisation& p ) const {
    const util::Config* params = dynamic_cast<const util::Config*>( &p );
    if ( !params ) {
        throw_Exception( "Parametrisation has to be of atlas::util::Config type" );
    }

    util::Config geometry;
    if ( !params->get( "geometry", geometry ) ) {
        throw_Exception( "Could not find 'geometry' in Parametrisation", Here() );
    }

    std::string grid_uid;
    if ( geometry.get( "grid", grid_uid ) ) {
        Grid grid( grid_uid );
        if ( !geometry.has( "ngptot" ) ) {
            geometry.set( "ngptot", grid.size() );
        }
    }

    if ( !geometry.has( "ngptot" ) ) {
        throw_Exception( "Could not find 'ngptot' in Parametrisation" );
    }

    std::vector<util::Config> fields;
    if ( params->get( "fields", fields ) ) {
        for ( size_t i = 0; i < fields.size(); ++i ) {
            util::Config fieldparams;
            // Subsequent "set" calls can overwrite eachother, so that
            // finetuning is possible in e.g. the fields Parametrisation (such as
            // creator, nlev, ngptot)
            fieldparams.set( "creator", "IFS" );
            fieldparams.set( geometry );
            fieldparams.set( fields[i] );
            state.add( Field( fieldparams ) );

            // debug info
            std::stringstream s;
            eckit::JSON json( s );
            json << fieldparams;
            Log::debug() << "fieldparams = " << s.str() << std::endl;
        }
    }
}

// Register in factory
StateGeneratorBuilder<MyStateGenerator> __MyStateGenerator( "MyStateGenerator" );

// ===================================================================
//                               BEGIN TESTS
// ===================================================================

CASE( "state" ) {
    State state;
    EXPECT( state.size() == 0 );

    state.add( Field( "myfield", array::make_datatype<double>(), array::make_shape( 10, 1 ) ) );
    state.add( Field( "", array::make_datatype<double>(), array::make_shape( 10, 2 ) ) );
    state.add( Field( "", array::make_datatype<double>(), array::make_shape( 10, 3 ) ) );

    EXPECT( state.size() == 3 );
    EXPECT( state.has( "myfield" ) );
    EXPECT( state.has( "field_00001" ) );
    EXPECT( state.has( "field_00002" ) );

    EXPECT( state.field( 0 ).name() == std::string( "field_00001" ) );
    EXPECT( state.field( 1 ).name() == std::string( "field_00002" ) );
    EXPECT( state.field( 2 ).name() == std::string( "myfield" ) );

    state.remove( "myfield" );
    EXPECT( state.size() == 2 );
    EXPECT( !state.has( "myfield" ) );

    state.remove( "field_00002" );
    EXPECT( state.size() == 1 );
    EXPECT( !state.has( "field_00002" ) );
}

CASE( "state_generator" ) {
    EXPECT( StateGeneratorFactory::has( "MyStateGenerator" ) );
    std::unique_ptr<StateGenerator> stategenerator( StateGeneratorFactory::build( "MyStateGenerator" ) );
}

CASE( "state_create" ) {
    util::Config p;
    util::Config geometry;
    geometry.set( "grid", "O80" );
    geometry.set( "ngptot", 350 );
    geometry.set( "nproma", 3 );
    geometry.set( "nlev", 5 );
    p.set( "geometry", geometry );

    std::vector<util::Config> fields( 5 );
    fields[0].set( "name", "temperature" );
    fields[0].set( "datatype", array::DataType::real32().str() );

    fields[1].set( "name", "wind" );
    fields[1].set( "nvar", 2 );  // vector field u,v
    fields[1].set( "datatype", array::DataType::real64().str() );

    fields[2].set( "name", "soiltype" );
    fields[2].set( "datatype", array::DataType::int32().str() );
    fields[2].set( "nlev", 1 );  // We can overwrite nlev from geometry here

    fields[3].set( "name", "GFL" );
    fields[3].set( "nvar", 12 );  // assume 12 variables in GFL array
    fields[3].set( "datatype", array::DataType::real64().str() );

    fields[4].set( "name", "array" );
    fields[4].set( "datatype", array::DataType::int64().str() );
    fields[4].set( "creator", "ArraySpec" );
    fields[4].set( "shape", array::make_shape( 10, 2 ) );

    p.set( "fields", fields );

    // We can also translate parameters to a json:
    std::stringstream json;
    eckit::JSON js( json );
    js << p;
    Log::info() << "json = " << json.str() << std::endl;

    // And we can create back parameters from json:
    util::Config from_json_stream( json );

    // And if we have a json file, we could create Parameters from the file:
    // StateGenerater::Parameters from_json_file( eckit::PathName("file.json") );

    State state( "MyStateGenerator", p );

    EXPECT( state.has( "temperature" ) );
    EXPECT( state.has( "wind" ) );
    EXPECT( state.has( "soiltype" ) );

    Log::info() << state.field( "temperature" ) << std::endl;
    Log::info() << state.field( "wind" ) << std::endl;
    Log::info() << state.field( "soiltype" ) << std::endl;
    Log::info() << state.field( "GFL" ) << std::endl;

    array::ArrayView<float, 4> temperature = array::make_view<float, 4>( state.field( "temperature" ) );
    temperature( 0, 0, 0, 0 )              = 0;

    array::ArrayView<double, 4> wind = array::make_view<double, 4>( state.field( "wind" ) );
    wind( 0, 0, 0, 0 )               = 0;

    array::ArrayView<int, 4> soiltype = array::make_view<int, 4>( state.field( "soiltype" ) );
    soiltype( 0, 0, 0, 0 )            = 0;

    array::ArrayView<long, 2> array = array::make_view<long, 2>( state["array"] );
    array( 0, 0 )                   = 0;
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
