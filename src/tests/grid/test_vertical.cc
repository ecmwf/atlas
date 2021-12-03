/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <vector>
#include "atlas/grid/Vertical.h"
#include "atlas/option.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

std::vector<double> zrange(idx_t nlev, double min, double max) {
    std::vector<double> zcoord(nlev);
    double dzcoord = (max - min) / double(nlev - 1);
    for (idx_t jlev = 0; jlev < nlev; ++jlev) {
        zcoord[jlev] = min + jlev * dzcoord;
    }
    return zcoord;
}

//-----------------------------------------------------------------------------

CASE("test vertical;  default constructor") {
    Vertical vertical;
    EXPECT(vertical.size() == 0);
    EXPECT(vertical.k_begin() == 0);
    EXPECT(vertical.k_end() == 0);
}
CASE("test vertical;  config levels") {
    Vertical vertical(option::levels(10));
    EXPECT(vertical.size() == 10);
    EXPECT(vertical.k_begin() == 0);
    EXPECT(vertical.k_end() == 10);
}
CASE("test vertical;  array") {
    Vertical vertical(5, zrange(5, 0., 1.));
    EXPECT(vertical.size() == 5);
    EXPECT(vertical.k_begin() == 0);
    EXPECT(vertical.k_end() == 5);
}

//-----------------------------------------------------------------------------


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
