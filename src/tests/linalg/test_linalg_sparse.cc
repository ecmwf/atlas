/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "tests/AtlasTestEnvironment.h"
#include "tests/linalg/helper_linalg_sparse.h"


using namespace atlas::linalg;

namespace atlas {
namespace test {

using SparseMatrix = eckit::linalg::SparseMatrix;

//----------------------------------------------------------------------------------------------------------------------

// strings to be used in the tests
static std::string eckit_linalg = sparse::backend::eckit_linalg::type();
static std::string openmp       = sparse::backend::openmp::type();
static std::string hicsparse    = sparse::backend::hicsparse::type();

//----------------------------------------------------------------------------------------------------------------------

CASE("test introspection") {
    SECTION("ArrayView") {
        array::ArrayT<float> array(4, 3);
        auto view  = array::make_view<const float, 2>(array);
        using View = decltype(view);
        static_assert(introspection::has_contiguous<View>::value, "ArrayView expected to have contiguous");
        static_assert(introspection::has_rank<View>::value, "ArrayView expected to have rank");
        static_assert(introspection::has_shape<View>::value, "ArrayView expected to have shape");
        EXPECT(introspection::contiguous(view));
        EXPECT_EQ(introspection::shape<0>(view), 4);
        EXPECT_EQ(introspection::shape<1>(view), 3);
    }
    SECTION("std::vector") {
        using View = std::vector<double>;
        static_assert(not introspection::has_contiguous<View>::value, "std::vector does not have contiguous");
        static_assert(not introspection::has_rank<View>::value, "std::vector does not have rank");
        static_assert(not introspection::has_shape<View>::value, "std::vector does not have shape");
        static_assert(introspection::rank<View>() == 1, "std::vector is of rank 1");
        auto v = View{1, 2, 3, 4};
        EXPECT(introspection::contiguous(v));
        EXPECT_EQ(introspection::shape<0>(v), 4);
        EXPECT_EQ(introspection::shape<1>(v), 4);
    }
    SECTION("eckit::linlag::Vector") {
        using View = Vector;
        static_assert(not introspection::has_contiguous<View>::value, "eckit::linalg::Vector does not have contiguous");
        static_assert(not introspection::has_rank<View>::value, "eckit::linalg::Vector does not have rank");
        static_assert(not introspection::has_shape<View>::value, "seckit::linalg::Vector does not have shape");
        static_assert(introspection::rank<View>() == 1, "eckit::linalg::Vector is of rank 1");
        auto v = Vector{1, 2, 3, 4};
        EXPECT(introspection::contiguous(v));
        EXPECT_EQ(introspection::shape<0>(v), 4);
    }
    SECTION("eckit::linlag::Matrix") {
        using View = Matrix;
        static_assert(not introspection::has_contiguous<View>::value, "eckit::linalg::Matrix does not have contiguous");
        static_assert(not introspection::has_rank<View>::value, "eckit::linalg::Matrix does not have rank");
        static_assert(not introspection::has_shape<View>::value, "seckit::linalg::Matrix does not have shape");
        static_assert(introspection::rank<View>() == 2, "eckit::linalg::Matrix is of rank 1");
        auto m = Matrix{{1., 2.}, {3., 4.}, {5., 6.}};
        EXPECT(introspection::contiguous(m));
        // Following is reversed because M has column-major ordering
        EXPECT_EQ(introspection::shape<0>(m), 3);
        EXPECT_EQ(introspection::shape<1>(m), 2);
    }
}

//----------------------------------------------------------------------------------------------------------------------

CASE("test backend functionalities") {
    sparse::current_backend(openmp);
    EXPECT_EQ(sparse::current_backend().type(), openmp);
    EXPECT_EQ(sparse::current_backend().getString("backend", "undefined"), "undefined");

    sparse::current_backend(eckit_linalg);
    EXPECT_EQ(sparse::current_backend().type(), "eckit_linalg");
    EXPECT_EQ(sparse::current_backend().getString("backend", "undefined"), "undefined");
    sparse::current_backend().set("backend", "default");
    EXPECT_EQ(sparse::current_backend().getString("backend"), "default");

    sparse::current_backend(openmp);
    EXPECT_EQ(sparse::current_backend().getString("backend", "undefined"), "undefined");
    EXPECT_EQ(sparse::default_backend(eckit_linalg).getString("backend"), "default");

    sparse::current_backend(hicsparse);
    EXPECT_EQ(sparse::current_backend().type(), "hicsparse");
    EXPECT_EQ(sparse::current_backend().getString("backend", "undefined"), "undefined");

    sparse::default_backend(eckit_linalg).set("backend", "generic");
    EXPECT_EQ(sparse::default_backend(eckit_linalg).getString("backend"), "generic");

    const sparse::Backend backend_default      = sparse::Backend();
    const sparse::Backend backend_openmp       = sparse::backend::openmp();
    const sparse::Backend backend_eckit_linalg = sparse::backend::eckit_linalg();
    const sparse::Backend backend_hicsparse    = sparse::backend::hicsparse();
    EXPECT_EQ(backend_default.type(), hicsparse);
    EXPECT_EQ(backend_openmp.type(), openmp);
    EXPECT_EQ(backend_eckit_linalg.type(), eckit_linalg);
    EXPECT_EQ(backend_hicsparse.type(), hicsparse);

    EXPECT_EQ(std::string(backend_openmp), openmp);
    EXPECT_EQ(std::string(backend_eckit_linalg), eckit_linalg);
}

//----------------------------------------------------------------------------------------------------------------------

CASE("sparse_matrix vector multiply (spmv)") {
    // "square" matrix
    // A =  2  . -3
    //      .  2  .
    //      .  .  2
    // x = 1 2 3
    // y = 1 2 3
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};

    for (std::string backend : {openmp, eckit_linalg}) {
        sparse::current_backend(backend);

        SECTION("test_identity [backend=" + sparse::current_backend().type() + "]") {
            {
                Vector y1(3);
                SparseMatrix B;
                B.setIdentity(3, 3);
                sparse_matrix_multiply(B, Vector{1., 2., 3.}, y1);
                expect_equal(y1, Vector{1., 2., 3.});
            }

            {
                SparseMatrix C;
                C.setIdentity(6, 3);
                Vector y2(6);
                sparse_matrix_multiply(C, Vector{1., 2., 3.}, y2);
                expect_equal(y2, Vector{1., 2., 3.});
                expect_equal(y2.data() + 3, Vector{0., 0., 0.}.data(), 3);
            }

            {
                SparseMatrix D;
                D.setIdentity(2, 3);
                Vector y3(2);
                sparse_matrix_multiply(D, Vector{1., 2., 3.}, y3);
                expect_equal(y3, Vector{1., 2., 3.});
            }
        }


        SECTION("eckit::Vector [backend=" + sparse::current_backend().type() + "]") {
            Vector y(3);
            sparse_matrix_multiply(A, Vector{1., 2., 3.}, y);
            expect_equal(y, Vector{-7., 4., 6.});
            // spmv of sparse matrix and vector of non-matching sizes should fail
            EXPECT_THROWS_AS(sparse_matrix_multiply(A, Vector(2), y), eckit::AssertionFailed);
        }

        SECTION("View of atlas::Array [backend=" + backend + "]") {
            ArrayVector<double> x(Vector{1., 2., 3.});
            ArrayVector<double> y(3);
            auto x_view = x.view();
            auto y_view = y.view();
            sparse_matrix_multiply(A, x_view, y_view);
            expect_equal(y.view(), Vector{-7., 4., 6.});
            // sparse_matrix_multiply of sparse matrix and vector of non-matching sizes should fail
            {
                ArrayVector<double> x2(2);
                auto x2_view = x2.view();
                EXPECT_THROWS_AS(sparse_matrix_multiply(A, x2_view, y_view), eckit::AssertionFailed);
            }
        }

        SECTION("View of atlas::Array [backend=" + backend + "]") {
            ArrayVector<double> x(Vector{1., 2., 3.});
            ArrayVector<double> y(3);
            auto x_view = x.view();
            auto y_view = y.view();
            auto spmm = SparseMatrixMultiply{backend};
            spmm(A, x_view, y_view);
            expect_equal(y_view, Vector{-7., 4., 6.});
        }

        SECTION("View of atlas::Array [backend=" + backend + "]") {
            ArrayVector<double> x(Vector{1., 2., 3.});
            ArrayVector<double> y(3);
            auto x_view = x.view();
            auto y_view = y.view();
            auto spmm = SparseMatrixMultiply{backend};
            spmm.multiply(A, x_view, y_view);
            expect_equal(y_view, Vector{-7., 4., 6.});
        }
    }
}

CASE("sparse_matrix vector multiply-add (spmv)") {
    // "square" matrix
    // A =  2  . -3
    //      .  2  .
    //      .  .  2
    // x = 1 2 3
    // y = 1 2 3
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};

    for (std::string backend : {openmp, eckit_linalg}) {
        sparse::current_backend(backend);

        SECTION("View of atlas::Array [backend=" + backend + "]") {
            ArrayVector<double> x(Vector{1., 2., 3.});
            ArrayVector<double> y(Vector{4., 5., 6.});
            auto x_view = x.view();
            auto y_view = y.view();
            sparse_matrix_multiply_add(A, x_view, y_view);
            expect_equal(y.view(), Vector{-3., 9., 12.});
            // sparse_matrix_multiply_add of sparse matrix and vector of non-matching sizes should fail
            {
                ArrayVector<double> x2(2);
                auto x2_view = x2.view();
                EXPECT_THROWS_AS(sparse_matrix_multiply_add(A, x2_view, y_view), eckit::AssertionFailed);
            }
        }

        SECTION("sparse_matrix_multiply_add [backend=" + backend + "]") {
            ArrayVector<double> x(Vector{1., 2., 3.});
            ArrayVector<double> y(Vector{1., 2., 3.});
            auto x_view = x.view();
            auto y_view = y.view();
            auto spmm = SparseMatrixMultiply{sparse::backend::openmp()};
            spmm.multiply_add(A, x_view, y_view);
            expect_equal(y_view, Vector{-6., 6., 9.});
        }
    }
}

CASE("sparse_matrix matrix multiply (spmm)") {
    // "square"
    // A =  2  . -3
    //      .  2  .
    //      .  .  2
    // x = 1 2 3
    // y = 1 2 3
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    Matrix m{{1., 2.}, {3., 4.}, {5., 6.}};
    Matrix c_exp{{-13., -14.}, {6., 8.}, {10., 12.}};

    for (std::string backend : {openmp, eckit_linalg}) {
        sparse::current_backend(backend);

        SECTION("eckit::Matrix [backend=" + sparse::current_backend().type() + "]") {
            auto c = Matrix(3, 2);
            sparse_matrix_multiply(A, m, c);
            expect_equal(c, c_exp);
            // spmm of sparse matrix and matrix of non-matching sizes should fail
            EXPECT_THROWS_AS(sparse_matrix_multiply(A, Matrix(2, 2), c), eckit::AssertionFailed);
        }


        SECTION("View of eckit::Matrix [backend=" + backend + "]") {
            auto c  = Matrix(3, 2);
            auto mv = atlas::linalg::make_view(m);  // convert Matrix to array::LocalView<double,2>
            auto cv = atlas::linalg::make_view(c);
            sparse_matrix_multiply(A, mv, cv, Indexing::layout_right);
            expect_equal(c, c_exp);
        }

        SECTION("View of atlas::Array PointsRight [backend=" + sparse::current_backend().type() + "]") {
            ArrayMatrix<double, Indexing::layout_right> ma(m);
            ArrayMatrix<double, Indexing::layout_right> c(3, 2);
            auto ma_view = ma.view();
            auto c_view  = c.view();
            sparse_matrix_multiply(A, ma_view, c_view, Indexing::layout_right);
            expect_equal(c_view, c_exp);
        }
    }

    SECTION("sparse_matrix_multiply [backend=openmp]") {
        sparse::current_backend(eckit_linalg);  // expected to be ignored further
        auto backend = sparse::backend::openmp();
        ArrayMatrix<float> ma(m);
        ArrayMatrix<float> c(3, 2);
        auto ma_view = ma.view();
        auto c_view  = c.view();
        sparse_matrix_multiply(A, ma_view, c_view, backend);
        expect_equal(c_view, ArrayMatrix<float>(c_exp).view());
    }

    SECTION("SparseMatrixMultiply [backend=openmp] 1") {
        sparse::current_backend(eckit_linalg);  // expected to be ignored
        auto spmm = SparseMatrixMultiply{sparse::backend::openmp()};
        ArrayMatrix<float> ma(m);
        ArrayMatrix<float> c(3, 2);
        auto ma_view = ma.view();
        auto c_view  = c.view();
        spmm(A, ma_view, c_view);
        expect_equal(c_view, ArrayMatrix<float>(c_exp).view());
    }

    SECTION("SparseMatrixMultiply [backend=openmp] 2") {
        sparse::current_backend(eckit_linalg);  // expected to be ignored
        auto spmm = SparseMatrixMultiply{openmp};
        ArrayMatrix<float> ma(m);
        ArrayMatrix<float> c(3, 2);
        auto ma_view = ma.view();
        auto c_view  = c.view();
        spmm(A, ma_view, c_view);
        expect_equal(c_view, ArrayMatrix<float>(c_exp).view());
    }

    SECTION("SparseMatrixMultiply::multiply [backend=openmp]") {
        sparse::current_backend(eckit_linalg);  // expected to be ignored
        auto spmm = SparseMatrixMultiply{openmp};
        ArrayMatrix<float> ma(m);
        ArrayMatrix<float> c(3, 2);
        auto ma_view = ma.view();
        auto c_view  = c.view();
        spmm.multiply(A, ma_view, c_view);
        expect_equal(c_view, ArrayMatrix<float>(c_exp).view());
    }
}

CASE("sparse_matrix matrix multiply-add (spmm)") {
    // "square"
    // A =  2  . -3
    //      .  2  .
    //      .  .  2
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    Matrix m{{1., 2.}, {3., 4.}, {5., 6.}};
    Matrix y_exp{{-12., -12.}, {9., 12.}, {15., 18.}};

    for (std::string backend : {openmp, eckit_linalg}) {
        sparse::current_backend(backend);

        SECTION("View of atlas::Array PointsRight [backend=" + sparse::current_backend().type() + "]") {
            ArrayMatrix<double, Indexing::layout_right> x(m);
            ArrayMatrix<double, Indexing::layout_right> y(m);
            auto x_view = x.view();
            auto y_view = y.view();
            sparse_matrix_multiply_add(A, x_view, y_view, Indexing::layout_right);
            expect_equal(y_view, y_exp);
        }
    }

    SECTION("sparse_matrix_multiply_add [backend=openmp]") {
        ArrayMatrix<double> x(m);
        ArrayMatrix<double> y(m);
        auto x_view = x.view();
        auto y_view = y.view();
        sparse_matrix_multiply_add(A, x_view, y_view, sparse::backend::openmp());
        expect_equal(y_view, ArrayMatrix<double>(y_exp).view());
    }

    SECTION("SparseMatrixMultiply::multiply_add [backend=openmp]") {
        auto spmm = SparseMatrixMultiply{sparse::backend::openmp()};
        ArrayMatrix<double> x(m);
        ArrayMatrix<double> y(m);
        auto x_view = x.view();
        auto y_view = y.view();
        spmm.multiply_add(A, x_view, y_view);
        expect_equal(y_view, ArrayMatrix<double>(y_exp).view());
    }
}


//----------------------------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
