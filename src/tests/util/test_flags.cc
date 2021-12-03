/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iostream>

#include "atlas/util/Bitflags.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::util::Bitflags;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_Flags") {
    int b1 = (1 << 0);
    int b2 = (1 << 1);
    int b3 = (1 << 2);
    int b4 = (1 << 3);

    int bits = b1 | b2;
    std::cout << Bitflags::bitstr(bits) << std::endl;
    EXPECT(bits == 3);

    EXPECT(Bitflags::check(bits, b1));
    EXPECT(Bitflags::check(bits, b2));
    EXPECT(!Bitflags::check(bits, b3));
    EXPECT(Bitflags::check_all(bits, b1 | b2));
    EXPECT(!Bitflags::check_all(bits, b1 | b3));
    EXPECT(Bitflags::check_any(bits, b1 | b3));
    EXPECT(!Bitflags::check_any(bits, b3 | b4));
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
