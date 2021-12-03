/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/MicroDeg.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("convert to microdegrees") {
    EXPECT_EQ(microdeg(-179.999999999999), -180000000);
    EXPECT_EQ(microdeg(+179.999999999999), +180000000);
    EXPECT_EQ(microdeg(-180.000000000001), -180000000);
    EXPECT_EQ(microdeg(+180.000000000001), +180000000);
    EXPECT_EQ(microdeg(-0.000000000001), 0);
    EXPECT_EQ(microdeg(+0.000000000001), 0);
    EXPECT_EQ(microdeg(-0.000001), -1);
    EXPECT_EQ(microdeg(+0.000001), +1);
    EXPECT_EQ(microdeg(-0.000001499999), -1);
    EXPECT_EQ(microdeg(+0.000001499999), +1);
    EXPECT_EQ(microdeg(-0.000001500001), -2);
    EXPECT_EQ(microdeg(+0.000001500001), +2);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
