/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"
#include "eckit/runtime/Tool.h"

#include "atlas/mpi/mpi.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Grid.h"
#include "atlas/grids/grids.h"
#include "atlas/FieldSet.h"
#include "atlas/meshgen/Tesselation.h"
#include "atlas/meshgen/Delaunay.h"
#include "atlas/io/Gmsh.h"

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
};

//-----------------------------------------------------------------------------

void TestField::test_constructor()
{
  // create a grid

  BoundBox earth ( Grid::Point(0.,-90.), Grid::Point(359.999999,90.) );

  Grid::Ptr g (new atlas::grids::LonLatGrid( 20, 10, earth ) );
//  Grid::Ptr g (Grid::create("oct.N6"));

  ASSERT( g );

  // Build a mesh from grid
  Mesh mesh(*g);

  ASSERT( mesh.grid().same( *g ) );
  ASSERT( mesh.has_function_space("nodes") );

  // create some reference data for testing

  std::vector<double> ref_data;
  ref_data.reserve( g->npts() );
  for(size_t i = 0; i < ref_data.size(); ++i)
    ref_data.push_back( (double)i );

  // now build a test field handle

  std::string sname("field_name");

  atlas::FunctionSpace& nodes  = mesh.function_space( "nodes" );

  FieldT<double>& data = nodes.create_field<double>( sname,1);

  for(size_t i = 0; i < ref_data.size(); i++)
    data[i] = ref_data[i];

  Field::Ptr f ( &nodes.field( sname ) );

  ASSERT( f );

  eckit::Log::info() << f.owners() << std::endl;

  ASSERT( f.owners() == 2 );

  Field::Vector fields;
  fields.push_back(f);

  atlas::FieldSet fs(fields);

  // iterate over the fields
  for (Field::Vector::iterator it = fs.fields().begin(); it != fs.fields().end(); ++it)
  {
    ArrayView<double> vdata( **it );

    for( size_t i = 0; i < ref_data.size(); ++i )
    {
      ASSERT( ref_data[i] == vdata(i) );
    }
  }
}

//-----------------------------------------------------------------------------

void TestField::run()
{
    eckit::mpi::init();
    test_constructor();
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

