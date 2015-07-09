/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include <string>

#include "atlas/atlas_config.h"

#define BOOST_TEST_MODULE test_state
#include "ecbuild/boost_test_framework.h"

#include "eckit/memory/ScopedPtr.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/parser/JSON.h"
#include "eckit/parser/JSONParser.h"

#include "atlas/atlas.h"
#include "atlas/State.h"
#include "atlas/Mesh.h"
#include "atlas/Field.h"
#include "atlas/Grid.h"
#include "atlas/util/DataType.h"
#include "atlas/util/ArrayView.h"

// ------------------------------------------------------------------

namespace atlas {
namespace test {

// -------------------------------------------------------------------
// Example State generator

// ---  Declaration (in .h file)
class MyStateGenerator : public StateGenerator {
public:
  MyStateGenerator( const eckit::Parametrisation& p = Parameters() ) : StateGenerator(p) {}
  ~MyStateGenerator() {}
  virtual void generate( State& state, const eckit::Parametrisation& p = Parameters() ) const;
};

// ---  Implementation (in .cc file)
void MyStateGenerator::generate( State& state, const eckit::Parametrisation& p ) const
{
  const atlas::Parametrisation *params = dynamic_cast<const atlas::Parametrisation*>(&p);
  if( !params )
  {
    throw eckit::Exception("Parametrisation has to be of atlas::Parametrisation type");
  }

  atlas::Parametrisation geometry;
  if( ! params->get("geometry",geometry) ) {
    throw eckit::BadParameter("Could not find 'geometry' in Parametrisation",Here());
  }
  
  std::string grid;
  if( geometry.get("grid",grid) )
  {
    state.add( Grid::create( grid ) );
    if( !geometry.has("ngptot") ) {
      geometry.set("ngptot",state.grid().npts());
    }
  }
  
  if( ! geometry.has("ngptot") ) {
    throw eckit::BadParameter("Could not find 'ngptot' in Parametrisation");
  }
  
  std::vector<atlas::Parametrisation> fields;
  if( params->get("fields",fields) )
  {
    for( size_t i=0; i<fields.size(); ++i )
    {
      atlas::Parametrisation fieldparams;
      // Subsequent "set" calls can overwrite eachother, so that
      // finetuning is possible in e.g. the fields Parametrisation (such as creator, nlev, ngptot)
      fieldparams.set("creator","IFS");
      fieldparams.set(geometry);
      fieldparams.set(fields[i]);
      state.add( Field::create( fieldparams ) );
    }    
  }
    //
  // eckit::ValueList fields = properties.get("fields");
  //
}

// Register in factory
StateGeneratorBuilder<MyStateGenerator> __MyStateGenerator("MyStateGenerator");

// ===================================================================
//                               BEGIN TESTS 
// ===================================================================


struct GlobalFixture {
    GlobalFixture()  { atlas_init(); }
    ~GlobalFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( GlobalFixture )

BOOST_AUTO_TEST_SUITE( test_state )

BOOST_AUTO_TEST_CASE( state )
{
  State state;
  BOOST_CHECK_EQUAL( state.nb_fields() , 0 );
  BOOST_CHECK_EQUAL( state.nb_meshes() , 0 );
  BOOST_CHECK_EQUAL( state.nb_grids()  , 0 );
  
  state.add( Field::create( make_shape(10,1) , Field::Parameters("name","myfield") ) );
  state.add( Field::create( make_shape(10,2) ) );
  state.add( Field::create( make_shape(10,3) ) );

  BOOST_CHECK_EQUAL( state.nb_fields() , 3 );
  BOOST_CHECK( state.has_field("myfield") );
  BOOST_CHECK( state.has_field("field_00001") );
  BOOST_CHECK( state.has_field("field_00002") );

  BOOST_CHECK_EQUAL( state.field(0).name(), std::string("field_00001") );
  BOOST_CHECK_EQUAL( state.field(1).name(), std::string("field_00002") );
  BOOST_CHECK_EQUAL( state.field(2).name(), std::string("myfield") );

  state.remove_field("myfield");
  BOOST_CHECK_EQUAL( state.nb_fields() , 2 );
  BOOST_CHECK( ! state.has_field("myfield") );

  state.remove_field("field_00002");
  BOOST_CHECK_EQUAL( state.nb_fields() , 1 );
  BOOST_CHECK( ! state.has_field("field_00002") );

  Grid& grid = state.add( Grid::create("rgg.N16") );
  Mesh& mesh = state.add( Mesh::create(grid) );
}

BOOST_AUTO_TEST_CASE( state_generator )
{
  BOOST_CHECK( StateGeneratorFactory::has("MyStateGenerator") );
  eckit::ScopedPtr<StateGenerator> stategenerator ( StateGeneratorFactory::build("MyStateGenerator") );
}

BOOST_AUTO_TEST_CASE( state_create )
{

  //-- Create JSON --------------
  eckit::Properties p;
  eckit::Properties geometry;
  geometry.set("grid","oct.N80");  
  geometry.set("ngptot",35000);
  geometry.set("nproma",20);  
  geometry.set("nlev",137);
  p.set("geometry",geometry);
  
  std::vector<eckit::Properties> fields(5);
  fields[0].set("name","temperature");
  fields[0].set("data_type",DataType::real32());
  
  fields[1].set("name","wind");
  fields[1].set("nvar",2); // vector field u,v
  fields[1].set("data_type",DataType::real64());
  
  fields[2].set("name","soiltype");
  fields[2].set("data_type",DataType::int32());
  fields[2].set("nlev",1); // We can overwrite nlev from geometry here

  fields[3].set("name","GFL");
  fields[3].set("nvar",30); // assume 30 variables in GFL array
  fields[3].set("data_type",DataType::real64());

  fields[4].set("name","array");
  fields[4].set("data_type",DataType::int64());
  fields[4].set("creator","ArraySpec");
  fields[4].set("shape",eckit::makeVectorValue( make_shape(10,2) ) );
  
  p.set("fields",eckit::makeVectorValue(fields));

  BOOST_CHECKPOINT("");
  
  std::stringstream json;
  eckit::JSON js(json);
  js << p;
  eckit::Log::info() << "JSON = " << json.str() << std::endl;
  //-- End Create JSON --------------
    
  eckit::ScopedPtr<State> state ( State::create("MyStateGenerator",StateGenerator::Parameters(json)) );
  
  BOOST_CHECK( state->has_field("temperature") );
  BOOST_CHECK( state->has_field("wind") );
  BOOST_CHECK( state->has_field("soiltype") );
  BOOST_CHECK( state->has_grid() );
  
  eckit::Log::info() << state->field("temperature") << std::endl;
  eckit::Log::info() << state->field("wind")        << std::endl;
  eckit::Log::info() << state->field("soiltype")    << std::endl;
  eckit::Log::info() << state->field("GFL")         << std::endl;
  
  ArrayView<float,4> temperature( state->field("temperature") );
  temperature(0,0,0,0) = 0;

  ArrayView<double,4> wind( state->field("wind") );
  wind(0,0,0,0) = 0;

  ArrayView<int,4> soiltype( state->field("soiltype") );
  soiltype(0,0,0,0) = 0;

  ArrayView<long,2> array( state->field("array") );
  array(0,0) = 0;

}

BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace atlas
