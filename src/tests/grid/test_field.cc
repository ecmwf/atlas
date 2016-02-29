/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/util/runtime/Log.h"
#include "eckit/runtime/Tool.h"
#include "eckit/value/CompositeParams.h"

#include "atlas/util/parallel/mpi/mpi.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/grid/grids.h"
#include "atlas/field/FieldSet.h"
#include "atlas/mesh/generators/Delaunay.h"
#include "atlas/util/io/Gmsh.h"
#include "atlas/internals/Debug.h"
#include "atlas/field/State.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/util/DataType.h"

using namespace std;
using namespace eckit;
using namespace atlas;
using namespace atlas::mesh::generators;

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
    void test_wrap_rawdata();
};

//-----------------------------------------------------------------------------

void TestField::test_constructor()
{
  // create a grid

  BoundBox earth ( Grid::Point(0.,-90.), Grid::Point(359.999999,90.) );

  Grid::Ptr g (new atlas::grid::LonLatGrid( 20ul, 10ul, earth ) );
//  Grid::Ptr g (Grid::create("O6"));

  ASSERT( g );

  // Build a mesh from grid
  mesh::Mesh mesh(*g);

  ASSERT( mesh.grid().same( *g ) );

  // create some reference data for testing

  std::vector<double> ref_data;
  ref_data.reserve( g->npts() );
  for(size_t i = 0; i < ref_data.size(); ++i)
    ref_data.push_back( (double)i );

  // now build a test field handle

  std::string sname("field_name");

  mesh::Nodes& nodes  = mesh.nodes();

  array::ArrayView<double,1> data ( nodes.add( field::Field::create<double>( sname,array::make_shape(nodes.size(),1) ) ) );

  for(size_t i = 0; i < ref_data.size(); i++)
    data[i] = ref_data[i];

  field::Field::Ptr f ( &nodes.field( sname ) );

  ASSERT( f );

  DEBUG_VAR( f.owners() );

  ASSERT( f.owners() == 2 );

  atlas::field::FieldSet fs;
  fs.add(*f);

  // iterate over the fields
  for (field::FieldSet::const_iterator it = fs.cbegin(); it != fs.cend(); ++it)
  {
    array::ArrayView<double> vdata( *it );

    for( size_t i = 0; i < ref_data.size(); ++i )
    {
      ASSERT( ref_data[i] == vdata(i) );
    }
  }
}




void TestField::test_fieldcreator()
{
  field::Field::Ptr field ( field::Field::create( util::Config
                                      ("creator","ArraySpec")
                                      ("shape",array::make_shape(10,2))
                                      ("datatype",DataType::real32().str())
                                      ("name","myfield")
                                  ));

  ASSERT( field->datatype() == DataType::real32() );
  ASSERT( field->name() == "myfield" );

  Grid::Ptr g (Grid::create("O6"));

  field::Field::Ptr arr (field::Field::create( util::Config
                                   ("creator","ArraySpec")
                                   ("shape",array::make_shape(10,2))
                               ));
  ASSERT( arr->shape(0) == 10 );
  ASSERT( arr->shape(1) == 2 );
  ASSERT( arr->datatype() == DataType::real64() );


  util::Config ifs_parameters;
  ifs_parameters
      ("creator","IFS")
      ("nlev",137)
      ("nproma",10)
      ("ngptot",g->npts());

  field::Field::Ptr ifs (field::Field::create( util::Config
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

  Log::debug() << std::flush;
  Log::info() << std::flush;
}

void take_array(const array::Array& arr)
{
  ASSERT( arr.size() == 20 );
}

class TakeArray
{
public:
  TakeArray(const array::Array& arr)
  {
    ASSERT( arr.size() == 20 );
  }
};

void TestField::test_implicit_conversion()
{
  SharedPtr<field::Field> field( field::Field::create<double>("tmp",array::make_shape(10,2)) );
  const array::Array& const_array = *field;
  array::Array& array = *field;

  array::ArrayView<double,2> arrv(array);
  arrv(0,0) = 8.;

  const array::ArrayView<double,2> carrv(const_array);
  ASSERT( carrv(0,0) == 8. );

  const array::ArrayView<double,2> cfieldv(*field);
  ASSERT( cfieldv(0,0) == 8. );

  take_array(*field);
  TakeArray ta(*field);

  const field::Field& f = *field;
  TakeArray cta(f);
}


void TestField::test_wrap_rawdata()
{
  std::vector<double> rawdata(20,8.);
  SharedPtr<array::Array> array( array::Array::wrap(rawdata.data(),array::make_shape(10,2)) );
  SharedPtr<field::Field> field( field::Field::create("wrapped",array.get()) );
  
  ASSERT( array->owners() == 2 );
  const array::ArrayView<double,2> cfieldv(*field);
  ASSERT( cfieldv(9,1) == 8. );
}


//-----------------------------------------------------------------------------

void TestField::run()
{
    eckit::mpi::init();
    test_constructor();
    test_fieldcreator();
    test_implicit_conversion();
    test_wrap_rawdata();
    eckit::mpi::finalize();
}

//-----------------------------------------------------------------------------

} // namespace test
} // namespace atlas

//-----------------------------------------------------------------------------

int main(int argc,char **argv)
{
    atlas::test::TestField t(argc,argv);
    return t.start();
}
