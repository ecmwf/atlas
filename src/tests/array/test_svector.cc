/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array/SVector.h"
#include "atlas/library/config.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::array;

namespace atlas {
namespace test {

CASE("test_svector") {
    SVector<int> list_ints(2);

    list_ints[0] = 3;
    list_ints[1] = 4;

    EXPECT(list_ints[0] == 3);
    EXPECT(list_ints[1] == 4);

    EXPECT(list_ints.size() == 2);

    list_ints[0]++;
    list_ints[1]++;

    EXPECT(list_ints[0] == 4);
    EXPECT(list_ints[1] == 5);
}

CASE("test_svector_resize") {
    SVector<int> list_ints(2);

    list_ints[0] = 3;
    list_ints[1] = 4;

    EXPECT(list_ints[0] == 3);
    EXPECT(list_ints[1] == 4);

    EXPECT(list_ints.size() == 2);

    list_ints.resize(5);

    EXPECT(list_ints.size() == 5);

    list_ints[3] = 5;
    list_ints[4] = 6;

    list_ints[3]++;
    list_ints[4]++;

    EXPECT(list_ints[0] == 3);
    EXPECT(list_ints[1] == 4);
    EXPECT(list_ints[3] == 6);
    EXPECT(list_ints[4] == 7);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
