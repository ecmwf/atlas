/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array/Array.h"
#include "atlas/array/ArrayViewUtil.h"
#include "atlas/array_fwd.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::array;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_var_size") {
    std::unique_ptr<Array> ds{Array::create<double>(4ul, 5ul, 7ul, 9ul)};

    auto arrv = make_host_view<double, 4>(*ds);
    EXPECT(ds->size() == 1260);

    EXPECT(get_var_size<0>(arrv) == 315);
    EXPECT(get_var_size<1>(arrv) == 252);
    EXPECT(get_var_size<2>(arrv) == 180);
    EXPECT(get_var_size<3>(arrv) == 140);
}

CASE("test_get_parallel_dim") {
    std::unique_ptr<Array> ds{Array::create<double>(4ul, 5ul, 7ul, 9ul)};

    auto arrv = make_host_view<double, 4>(*ds);
    EXPECT(ds->size() == 1260);

    EXPECT(get_parallel_dim<FirstDim>(arrv) == 0);
    EXPECT(get_parallel_dim<LastDim>(arrv) == 3);
    EXPECT(get_parallel_dim<Dim<1>>(arrv) == 1);
    EXPECT(get_parallel_dim<Dim<2>>(arrv) == 2);
}

CASE("test mixed index types") {
    std::unique_ptr<Array> ds{Array::create<double>(4ul, 5ul, 7ul, 9ul)};

    auto arrv = make_host_view<double, 4>(*ds);

    arrv(long(0), size_t(0), int(0), unsigned(0)) = 0.;
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
