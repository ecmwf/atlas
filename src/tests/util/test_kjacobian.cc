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

    // Default constructor.
    Jacobian2 a{};

    EXPECT_EQ(a.rows(), 2);
    EXPECT_EQ(a.cols(), 2);


    // Initialiser list constructor.
    const Jacobian2 b = {{-1., 1.5},
                         {1., -1.}};

    EXPECT_EQ(b.rows(), 2);
    EXPECT_EQ(b.cols(), 2);
    EXPECT_EQ(b(0, 0), -1.);
    EXPECT_EQ(b(0, 1), 1.5);
    EXPECT_EQ(b(1, 0), 1.);
    EXPECT_EQ(b(1, 1), -1.);


    // Base class constructor (base class methods return base class objects).
    const Jacobian2 c = b.inverse();

    EXPECT_EQ(c.rows(), 2);
    EXPECT_EQ(c.cols(), 2);
    EXPECT_APPROX_EQ(c(0, 0), 2.);
    EXPECT_APPROX_EQ(c(0, 1), 3.);
    EXPECT_APPROX_EQ(c(1, 0), 2.);
    EXPECT_APPROX_EQ(c(1, 1), 2.);

    // Trigger exception on list constructor.
    do {
        try {
            const Jacobian3 d = {{-1., 1.5}, {1., -1.}};                                             }
        catch (...) {
            break;
        }
        throw eckit::testing::TestException("Jacobian list constuctor missed error.", Here());
    } while (false);


    // Jacobian-Kpoint multiplication.
    const Jacobian3 A = {{0., -1., 0.},
                         {1., 0., 0.},
                         {0., 0., 1.}};

    const Point3 x = {1., 0., 0.};
    const Point3 Ax = {0., 1., 0.};

    EXPECT_EQ(A.rows(), 3);
    EXPECT_EQ(A.cols(), 3);
    EXPECT_EQ(A * x, Ax);


    // BaseMatrix-KPoint multiplication.
    const Point3 y = A.inverse().inverse() * x;

    EXPECT_APPROX_EQ(Ax[0], y[0]);
    EXPECT_APPROX_EQ(Ax[1], y[1]);
    EXPECT_APPROX_EQ(Ax[2], y[2]);

}

} // namespace test
} // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
