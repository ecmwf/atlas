/*
 * (C) Copyright 1996-2012 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"
#include "eckit/runtime/Tool.h"

#include "eckit/grid/Grid.h"
#include "eckit/grid/LatLon.h"
#include "eckit/grid/Field.h"

using namespace std;
using namespace eckit;

//-----------------------------------------------------------------------------

namespace eckit_test {

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
    using namespace eckit::grid;

    BoundBox2D earth ( Point2D(-90.,0.), Point2D(90.,360.) );
    Grid* g = new LatLon( 4, 4, earth );
    ASSERT( g );

    Field::MetaData* meta = new Field::MetaData();
    ASSERT( meta );

    Field::Data*     data = new Field::Data();
    ASSERT( data );

    // create some reference data for testing
    Field::Data ref_data;
    for (unsigned int i = 0; i < 1000; i++)
        ref_data.push_back((double)i);

    // copy the ref_data to the data that will be put into the
    // Field object
    for (unsigned int i = 0; i < ref_data.size(); i++)
        data->push_back(ref_data[i]);

    Field* f = new Field(g,meta,data);
    ASSERT( f );

    Field::Vector fields;
    fields.push_back(f);

    FieldSet fs(fields);
    //ASSERT(fs.grid().get() == g); 
    
    // iterate over the fields
    for (Field::Vector::iterator it = fs.fields().begin(); it != fs.fields().end(); ++it)
    {
        // extract and test the data
        Field::Data& d = (*it)->data();
        for (unsigned int i = 0; i < ref_data.size(); i++)
        {
            ASSERT(ref_data[i] == d[i]);
        }   
    }


}

//-----------------------------------------------------------------------------

void TestField::run()
{
    test_constructor();
}

//-----------------------------------------------------------------------------

} // namespace eckit_test

//-----------------------------------------------------------------------------

int main(int argc,char **argv)
{
    eckit_test::TestField mytest(argc,argv);
    mytest.start();
    return 0;
}

