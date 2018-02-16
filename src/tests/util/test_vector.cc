/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array/Vector.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::array;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE( "test_vector" ) {
    Vector<int> vec( 3 );

    VectorView<int> vec_view = make_host_vector_view( vec );

    vec_view[0] = 3;
    vec_view[1] = -3;
    vec_view[2] = 1;

    EXPECT( vec_view.size() == 3 );
    EXPECT( vec_view[0] == 3 );
    EXPECT( vec_view[1] == -3 );
    EXPECT( vec_view[2] == 1 );

    vec.resize( 5 );

    // TODO invalidate preview views
    VectorView<int> vec_viewb = make_host_vector_view( vec );
    vec_viewb[3]              = 5;
    vec_viewb[4]              = 6;

    EXPECT( vec_viewb[0] == 3 );
    EXPECT( vec_viewb[1] == -3 );
    EXPECT( vec_viewb[2] == 1 );
    EXPECT( vec_viewb[3] == 5 );
    EXPECT( vec_viewb[4] == 6 );
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main( int argc, char** argv ) {
    return atlas::test::run( argc, argv );
}
