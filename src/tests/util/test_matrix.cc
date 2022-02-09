/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iostream>

#include "atlas/util/Matrix.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using namespace util;

CASE("test matrix operations") {

    SECTION("addition") {

        const auto A = Matrix<double, 2, 2>{
            {1., 2.},
            {3., 4.}
        };

        const auto B = Matrix<double, 2, 2>{
            {5., 6.},
            {7., 8.}
        };

        const auto APlusB = Matrix<double, 2, 2>{
            {6., 8.},
            {10., 12.}
        };

        EXPECT(A + B == APlusB);
        EXPECT(APlusB - B == A);

        Log::info() << "A = " << std::endl << A << std::endl;
        Log::info() << "B = " << std::endl << B << std::endl;
        Log::info() << "A + B = " << std::endl << A + B << std::endl << std::endl;

    }

    SECTION("multiplication") {

        const auto A = Matrix<double, 4, 3>{
            {1., 0., 1.},
            {2., 1., 1.},
            {0., 1., 1.},
            {1., 1., 2.}
        };

        const auto B = Matrix<double, 3, 3>{
            {1., 2., 1.},
            {2., 3., 1.},
            {4., 2., 2.}
        };

        const auto AB = Matrix<double, 4, 3> {
            {5., 4., 3.},
            {8., 9., 5.},
            {6., 5., 3.},
            {11., 9., 6.}
        };

        EXPECT(A * B == AB);

        Log::info() << "A = " << std::endl << A << std::endl;
        Log::info() << "B = " << std::endl << B << std::endl;
        Log::info() << "A x B = " << std::endl << A * B << std::endl << std::endl;

    }

    SECTION("inverse") {

        const auto A = Matrix<double, 2, 2>{
            {-1., 1.5},
            {1., -1.}
        };

        const auto invA = Matrix<double, 2, 2>{
            {2., 3.},
            {2., 2.}
        };

        for (size_t i =0; i < A.size(); ++i) {
            EXPECT_APPROX_EQ(A.inverse().data()[i], invA.data()[i]);
        }

        Log::info() << "A = " << std::endl << A << std::endl;
        Log::info() << "A^-1 = " << std::endl << A.inverse() << std::endl << std::endl;

    }

    SECTION("transpose") {

        // Issue: transpose() is non const in eckit/maths/MatrixLapack.h!
        auto A = Matrix<double, 2, 2>{
            {1., 2.},
            {3., 4.}
        };


        const auto At = Matrix<double, 2, 2>{
            {1., 3.},
            {2., 4.}
        };

        EXPECT(A.transpose() == At);

        Log::info() << "A = " << std::endl << A << std::endl;
        Log::info() << "A^T = " << std::endl << A.transpose() << std::endl << std::endl;

    }

    SECTION("determinant") {

        const auto A = Matrix<double, 2, 2>{
            {3., 7.},
            {1., -4.}
        };

        EXPECT_APPROX_EQ(A.determinant(), -19.);

        Log::info() << "A = " << std::endl << A << std::endl;
        Log::info() << "|A| = " << std::endl << A.determinant() << std::endl << std::endl;

    }

    SECTION("Row/Col slice") {

        const auto A = Matrix<double, 3, 3>{
            {1., 2., 1.},
            {2., 3., 1.},
            {4., 2., 2.}
        };

        const auto ARow0 = Matrix<double, 1, 3>{
            {1., 2., 1.}
        };

        const auto ACol0 = Matrix<double, 3, 1>{
            {1.},
            {2.},
            {4.}
        };

        EXPECT(A.row(0) == ARow0);
        EXPECT(A.col(0) == ACol0);

        Log::info() << "A = " << std::endl << A << std::endl;
        Log::info() << "A(0, j) = " << std::endl << ARow0 << std::endl;
        Log::info() << "A(i, 0) = " << std::endl << ACol0 << std::endl << std::endl;

    }

    SECTION("Point linear transform") {

        const auto x = Point3{1., 0., 0.};

        const auto A = Matrix<double, 3, 3>{
            {0., -1., 0.},
            {1., 0., 0.},
            {0., 0., 1.}
        };

        const auto Ax = Point3{0., 1., 0.};

        EXPECT_EQ(A * x, Ax);

        Log::info() << "A = " << std::endl << A << std::endl;
        Log::info() << "x = " << std::endl << x << std::endl;
        Log::info() << "Ax = " << std::endl << Ax << std::endl << std::endl;

    }

    SECTION("Sign") {

        const auto A = Matrix<double, 2, 2>{
            {3., 7.},
            {0., -4.}
        };

        const auto signA = Matrix<double, 2, 2>{
            {1., 1.},
            {0., -1.}
        };

        EXPECT(A.sign() == signA);
        Log::info() << "A = " << std::endl << A << std::endl;
        Log::info() << "sign(A) = " << std::endl << signA << std::endl << std::endl;

    }

}

} // namespace test
} // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
