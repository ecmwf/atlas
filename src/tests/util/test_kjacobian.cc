/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>

#include "atlas/util/KJacobian.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using namespace util;

CASE("Test Jacobian") {

    // Make Jacobian object.
    const auto a = make_Jacobian2({{-1., 1.5},
                                   {1., -1.}});

    EXPECT_EQ(a.rows(), 2);
    EXPECT_EQ(a.cols(), 2);
    EXPECT_EQ(a(0, 0), -1.);
    EXPECT_EQ(a(0, 1), 1.5);
    EXPECT_EQ(a(1, 0), 1.);
    EXPECT_EQ(a(1, 1), -1.);

    // Trigger exception on make_Jacobian.
    do {
        try {
            const auto b = make_Jacobian3({{-1., 1.5},
                                           {1., -1.}});
        }
        catch (...) {
            break;
        }
        throw eckit::testing::TestException("Jacobian list constuctor missed error.", Here());
    } while (false);


    // Jacobian-Kpoint multiplication.
    const auto A = make_Jacobian3({{0., -1., 0.},
                                   {1., 0., 0.},
                                   {0., 0., 1.}});

    const Point3 x = {1., 0., 0.};
    const PointXYZ xyz = x;
    const Point3 Ax = {0., 1., 0.};

    EXPECT_EQ(A.rows(), 3);
    EXPECT_EQ(A.cols(), 3);
    EXPECT_EQ(A * x, Ax);
    EXPECT_EQ(A * xyz, Ax);

}

} // namespace test
} // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
