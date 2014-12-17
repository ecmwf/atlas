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

#include "atlas/mpl/MPL.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Grid.h"
#include "atlas/grids/LonLatGrid.h"
#include "atlas/FieldSet.h"
#include "atlas/Tesselation.h"

using namespace std;
using namespace eckit;
using namespace atlas;

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

	Grid::BoundBox earth ( Grid::Point(0.,-90.), Grid::Point(359.999999,90.) );

	Grid::Ptr g (new atlas::grids::LonLatGrid( 20, 10, earth ) );

    // create some reference data for testing

    std::vector<double> ref_data;
    ref_data.reserve( g->npts() );
    for(size_t i = 0; i < ref_data.size(); ++i)
		ref_data.push_back( (double)i );

    // now build a test field handle

    std::string sname("field_name");

	ASSERT( g );

	Mesh& mesh = g->mesh();

	ASSERT( mesh.grid().uid() == g->uid() );

    ASSERT( mesh.has_function_space("nodes") );

    atlas::FunctionSpace& nodes  = mesh.function_space( "nodes" );

    FieldT<double>& data = nodes.create_field<double>( sname,1);

	for(size_t i = 0; i < ref_data.size(); i++)
        data[i] = ref_data[i];

    // create field handle

	Field::Ptr f ( &nodes.field( sname ) );

	ASSERT( f );

	std::cout << f.owners() << std::endl;

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
    MPL::init();
    test_constructor();
    MPL::finalize();
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

