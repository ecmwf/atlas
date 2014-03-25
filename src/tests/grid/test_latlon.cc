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

using namespace std;
using namespace eckit;

//-----------------------------------------------------------------------------

namespace eckit_test {

//-----------------------------------------------------------------------------

class TestLatLon : public Tool {
public:

    TestLatLon(int argc,char **argv): Tool(argc,argv) {}

    ~TestLatLon() {}

    virtual void run();

    void test_constructor();
};

//-----------------------------------------------------------------------------

void TestLatLon::test_constructor()
{
    using namespace eckit::grid;

    BoundBox2D earth ( Point2D(-90.,0.), Point2D(90.,360.) );
    Grid* g = NULL;

    // standard case

    g = new LatLon( 4, 4, earth );

    ASSERT( g );
    ASSERT( g->coordinates().size() == 25 );

    /// @todo substitute these comparisons with proper floating point comparisons
    ASSERT( g->boundingBox().bottom_left_.lat_ == -90. );
    ASSERT( g->boundingBox().bottom_left_.lon_ ==   0. );
    ASSERT( g->boundingBox().top_right_.lat_ ==  90. );
    ASSERT( g->boundingBox().top_right_.lon_ == 360. );

    delete g; g = NULL;

    // 1x1 case

    g = new LatLon( 1, 1, earth );

    ASSERT( g );
    ASSERT( g->coordinates().size() == 4 );

    /// @todo substitute these comparisons with proper floating point comparisons
    ASSERT( g->boundingBox().bottom_left_.lat_ == -90. );
    ASSERT( g->boundingBox().bottom_left_.lon_ ==   0. );
    ASSERT( g->boundingBox().top_right_.lat_ ==  90. );
    ASSERT( g->boundingBox().top_right_.lon_ == 360. );

    delete g; g = NULL;

}

//-----------------------------------------------------------------------------

void TestLatLon::run()
{
    test_constructor();
}

//-----------------------------------------------------------------------------

} // namespace eckit_test

//-----------------------------------------------------------------------------

int main(int argc,char **argv)
{
    eckit_test::TestLatLon mytest(argc,argv);
    mytest.start();
    return 0;
}

