/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/array/Array.h"
#include "atlas/array/ArrayViewUtil.h"
#include "tests/AtlasTestEnvironment.h"
#include "eckit/testing/Test.h"

using namespace atlas::array;
using namespace eckit::testing;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_var_size") {

    Array* ds = Array::create<double>(4ul, 5ul, 7ul, 9ul);

    auto arrv = make_host_view<double, 4>(*ds);
    EXPECT(ds->size() == 1260);

    EXPECT( get_var_size<0>(arrv) == 315);
    EXPECT( get_var_size<1>(arrv) == 252);
    EXPECT( get_var_size<2>(arrv) == 180);
    EXPECT( get_var_size<3>(arrv) == 140);

    delete ds;

}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas


int main(int argc, char **argv) {
    atlas::test::AtlasTestEnvironment env( argc, argv );
    return run_tests ( argc, argv, false );
}

