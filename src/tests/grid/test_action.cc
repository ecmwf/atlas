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
#include "eckit/grid/Action.h"
#include "eckit/types/Types.h"

using namespace eckit;

//-----------------------------------------------------------------------------

namespace eckit_test {

//-----------------------------------------------------------------------------

class TestAction : public Tool {
public:

    TestAction(int argc,char **argv): Tool(argc,argv) {}

    ~TestAction() {}

    virtual void run();

    void test_execute();
};

//-----------------------------------------------------------------------------

void TestAction::test_execute()
{
    using namespace eckit::grid;

    BoundBox2D earth ( Point2D(-90.,0.), Point2D(90.,360.) );
    Grid* g_in = new LatLon( 4, 4, earth );
    ASSERT( g_in );

    FieldSet* fs_in = new FieldSet(g_in);
    
    Grid* g_out = new LatLon( 3, 3, earth );
    ASSERT( g_out );
    
    FieldSet* fs_out = new FieldSet(g_out);
    
    StringDict config;
    config["interpolator"] = "bilinear";
    Interpolate interpolate(config);

    interpolate( *fs_in, *fs_out );

    delete fs_in;
    delete fs_out;

}

//-----------------------------------------------------------------------------

void TestAction::run()
{
    test_execute();
}

//-----------------------------------------------------------------------------

} // namespace eckit_test

//-----------------------------------------------------------------------------

int main(int argc,char **argv)
{
    eckit_test::TestAction mytest(argc,argv);
    mytest.start();
    return 0;
}

