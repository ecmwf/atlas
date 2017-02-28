/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include <string>

#include "atlas/internals/atlas_config.h"

#define BOOST_TEST_MODULE test_state
#include "ecbuild/boost_test_framework.h"

#include "eckit/memory/ScopedPtr.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/parser/JSON.h"
#include "eckit/parser/JSONParser.h"

#include "atlas/atlas.h"
#include "atlas/field/State.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Grid.h"
#include "atlas/array/DataType.h"
#include "atlas/array/ArrayView.h"
#include "atlas/runtime/Log.h"

#include "tests/AtlasFixture.h"

using namespace atlas::field;

namespace atlas {
namespace test {

// -------------------------------------------------------------------
// Example State generator

// ---  Declaration (in .h file)
class MyStateGenerator : public StateGenerator {
public:
  MyStateGenerator( const eckit::Parametrisation& p = util::Config() ) : StateGenerator(p) {}
  ~MyStateGenerator() {}
  virtual void generate( State& state, const eckit::Parametrisation& p = util::Config() ) const;
};

// ---  Implementation (in .cc file)
void MyStateGenerator::generate( State& state, const eckit::Parametrisation& p ) const
{
  const util::Config *params = dynamic_cast<const util::Config*>(&p);
  if( !params )
  {
    throw eckit::Exception("Parametrisation has to be of atlas::util::Config type");
  }

  util::Config geometry;
  if( ! params->get("geometry",geometry) ) {
    throw eckit::BadParameter("Could not find 'geometry' in Parametrisation",Here());
  }

  std::string grid_uid;
  if( geometry.get("grid", grid_uid) )
  {
    grid::Grid grid(grid_uid);
    if (!geometry.has("ngptot")) {
      geometry.set("ngptot", grid.npts());
    }
  }

  if( ! geometry.has("ngptot") ) {
    throw eckit::BadParameter("Could not find 'ngptot' in Parametrisation");
  }

  std::vector<util::Config> fields;
  if( params->get("fields",fields) )
  {
    for( size_t i=0; i<fields.size(); ++i )
    {
      util::Config fieldparams;
      // Subsequent "set" calls can overwrite eachother, so that
      // finetuning is possible in e.g. the fields Parametrisation (such as creator, nlev, ngptot)
      fieldparams.set("creator","IFS");
      fieldparams.set(geometry);
      fieldparams.set(fields[i]);
      state.add( field::Field::create( fieldparams ) );

      // debug info
      std::stringstream s;
      eckit::JSON json(s);
      json << fieldparams;
      Log::debug() << "fieldparams = " << s.str() << std::endl;
    }
  }
}

// Register in factory
StateGeneratorBuilder<MyStateGenerator> __MyStateGenerator("MyStateGenerator");

// ===================================================================
//                               BEGIN TESTS
// ===================================================================

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_SUITE( test_state )

BOOST_AUTO_TEST_CASE( state )
{
  State state;
  BOOST_CHECK_EQUAL( state.size() , 0 );

  state.add( field::Field::create<double>( "myfield", array::make_shape(10,1) ) );
  state.add( field::Field::create<double>( "", array::make_shape(10,2) ) );
  state.add( field::Field::create<double>( "", array::make_shape(10,3) ) );

  BOOST_CHECK_EQUAL( state.size() , 3 );
  BOOST_CHECK( state.has("myfield") );
  BOOST_CHECK( state.has("field_00001") );
  BOOST_CHECK( state.has("field_00002") );

  BOOST_CHECK_EQUAL( state.field(0).name(), std::string("field_00001") );
  BOOST_CHECK_EQUAL( state.field(1).name(), std::string("field_00002") );
  BOOST_CHECK_EQUAL( state.field(2).name(), std::string("myfield") );

  state.remove("myfield");
  BOOST_CHECK_EQUAL( state.size() , 2 );
  BOOST_CHECK( ! state.has("myfield") );

  state.remove("field_00002");
  BOOST_CHECK_EQUAL( state.size() , 1 );
  BOOST_CHECK( ! state.has("field_00002") );

}

BOOST_AUTO_TEST_CASE( state_generator )
{
  BOOST_CHECK( StateGeneratorFactory::has("MyStateGenerator") );
  eckit::ScopedPtr<StateGenerator> stategenerator ( StateGeneratorFactory::build("MyStateGenerator") );
}

BOOST_AUTO_TEST_CASE( state_create )
{

  util::Config p;
  util::Config geometry;
  geometry.set("grid","O80");
  geometry.set("ngptot",350);
  geometry.set("nproma",3);
  geometry.set("nlev",5);
  p.set("geometry",geometry);

  std::vector<util::Config> fields(5);
  fields[0].set("name","temperature");
  fields[0].set("datatype",array::DataType::real32().str());

  fields[1].set("name","wind");
  fields[1].set("nvar",2); // vector field u,v
  fields[1].set("datatype",array::DataType::real64().str());

  fields[2].set("name","soiltype");
  fields[2].set("datatype",array::DataType::int32().str());
  fields[2].set("nlev",1); // We can overwrite nlev from geometry here

  fields[3].set("name","GFL");
  fields[3].set("nvar",12); // assume 12 variables in GFL array
  fields[3].set("datatype",array::DataType::real64().str());

  fields[4].set("name","array");
  fields[4].set("datatype",array::DataType::int64().str());
  fields[4].set("creator","ArraySpec");
  fields[4].set("shape",array::make_shape(10,2));

  p.set("fields",fields);

  // We can also translate parameters to a json:
  std::stringstream json;
  eckit::JSON js(json); js << p;
  Log::info() << "json = " << json.str() << std::endl;

  // And we can create back parameters from json:
  util::Config from_json_stream(json);

  // And if we have a json file, we could create Parameters from the file:
    // StateGenerater::Parameters from_json_file( eckit::PathName("file.json") );

  State state ( "MyStateGenerator",p );

  BOOST_CHECK( state.has("temperature") );
  BOOST_CHECK( state.has("wind") );
  BOOST_CHECK( state.has("soiltype") );

  Log::info() << state.field("temperature") << std::endl;
  Log::info() << state.field("wind")        << std::endl;
  Log::info() << state.field("soiltype")    << std::endl;
  Log::info() << state.field("GFL")         << std::endl;

  array::ArrayView<float,4> temperature( state.field("temperature") );
  temperature(0,0,0,0) = 0;

  array::ArrayView<double,4> wind( state.field("wind") );
  wind(0,0,0,0) = 0;

  array::ArrayView<int,4> soiltype( state.field("soiltype") );
  soiltype(0,0,0,0) = 0;

  array::ArrayView<long,2> array( state["array"] );
  array(0,0) = 0;

}

BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace atlas
