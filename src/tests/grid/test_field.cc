/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"
#include "eckit/runtime/Tool.h"
#include "eckit/value/CompositeParams.h"

#include "atlas/mpi/mpi.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Grid.h"
#include "atlas/Nodes.h"
#include "atlas/grids/grids.h"
#include "atlas/FieldSet.h"
#include "atlas/meshgen/Delaunay.h"
#include "atlas/io/Gmsh.h"
#include "atlas/util/Debug.h"
#include "atlas/State.h"
#include "atlas/Mesh.h"
#include "atlas/util/DataType.h"

using namespace std;
using namespace eckit;
using namespace atlas;
using namespace atlas::meshgen;

//-----------------------------------------------------------------------------

namespace atlas {
namespace  test {

//-----------------------------------------------------------------------------

class TestField : public Tool {
public:

    TestField(int argc,char **argv): Tool(argc,argv) {}

    ~TestField() {}

    virtual void run();

    void test_constructor();
    void test_fieldcreator();
    void test_implicit_conversion();
};

//-----------------------------------------------------------------------------

void TestField::test_constructor()
{
  // create a grid

  BoundBox earth ( Grid::Point(0.,-90.), Grid::Point(359.999999,90.) );

  Grid::Ptr g (new atlas::grids::LonLatGrid( 20ul, 10ul, earth ) );
//  Grid::Ptr g (Grid::create("O6"));

  ASSERT( g );

  // Build a mesh from grid
  Mesh mesh(*g);

  ASSERT( mesh.grid().same( *g ) );

  // create some reference data for testing

  std::vector<double> ref_data;
  ref_data.reserve( g->npts() );
  for(size_t i = 0; i < ref_data.size(); ++i)
    ref_data.push_back( (double)i );

  // now build a test field handle

  std::string sname("field_name");

  atlas::Nodes& nodes  = mesh.nodes();

  ArrayView<double,1> data ( nodes.add( Field::create<double>( sname,make_shape(nodes.size(),1) ) ) );

  for(size_t i = 0; i < ref_data.size(); i++)
    data[i] = ref_data[i];

  Field::Ptr f ( &nodes.field( sname ) );

  ASSERT( f );

  DEBUG_VAR( f.owners() );

  ASSERT( f.owners() == 2 );

  atlas::FieldSet fs;
  fs.add(*f);

  // iterate over the fields
  for (FieldSet::const_iterator it = fs.cbegin(); it != fs.cend(); ++it)
  {
    ArrayView<double> vdata( *it );

    for( size_t i = 0; i < ref_data.size(); ++i )
    {
      ASSERT( ref_data[i] == vdata(i) );
    }
  }
}




void TestField::test_fieldcreator()
{
  Field::Ptr field ( Field::create( Config
                                      ("creator","ArraySpec")
                                      ("shape",make_shape(10,2))
                                      ("datatype",DataType::real32().str())
                                      ("name","myfield")
                                  ));

  ASSERT( field->datatype() == DataType::real32() );
  ASSERT( field->name() == "myfield" );

  Grid::Ptr g (Grid::create("O6"));

  Field::Ptr arr (Field::create( Config
                                   ("creator","ArraySpec")
                                   ("shape",make_shape(10,2))
                               ));
  ASSERT( arr->shape(0) == 10 );
  ASSERT( arr->shape(1) == 2 );
  ASSERT( arr->datatype() == DataType::real64() );


  Config ifs_parameters;
  ifs_parameters
      ("creator","IFS")
      ("nlev",137)
      ("nproma",10)
      ("ngptot",g->npts());

  Field::Ptr ifs (Field::create( Config
                                    (ifs_parameters)
                                    ("name","myfield")
                                    ("datatype",DataType::int32().str())
                                    ("nvar",8)
                               ));

  DEBUG_VAR( *ifs );
  ASSERT( ifs->shape(0) == 36 );
  ASSERT( ifs->shape(1) == 8 );
  ASSERT( ifs->shape(2) == 137 );
  ASSERT( ifs->shape(3) == 10 );

  eckit::Log::debug() << std::flush;
  eckit::Log::info() << std::flush;
}

void take_array(const Array& arr)
{
  ASSERT( arr.size() == 20 );
}

class TakeArray
{
public:
  TakeArray(const Array& arr)
  {
    ASSERT( arr.size() == 20 );
  }
};

void TestField::test_implicit_conversion()
{
  SharedPtr<Field> field( Field::create<double>("tmp",make_shape(10,2)) );
  const Array& const_array = *field;
  Array& array = *field;

  ArrayView<double,2> arrv(array);
  arrv(0,0) = 8.;

  const ArrayView<double,2> carrv(const_array);
  ASSERT( carrv(0,0) == 8. );

  const ArrayView<double,2> cfieldv(*field);
  ASSERT( cfieldv(0,0) == 8. );

  take_array(*field);
  TakeArray ta(*field);

  const Field& f = *field;
  TakeArray cta(f);
}



//-----------------------------------------------------------------------------

void TestField::run()
{
    eckit::mpi::init();
    test_constructor();
    test_fieldcreator();
    test_implicit_conversion();
    eckit::mpi::finalize();
}

//-----------------------------------------------------------------------------

} // namespace test
} // namespace atlas

//-----------------------------------------------------------------------------

int main(int argc,char **argv)
{
	atlas::test::TestField mytest(argc,argv);
    mytest.start();
    return 0;
}

