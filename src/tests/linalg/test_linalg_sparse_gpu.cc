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
static std::string hicsparse = sparse::backend::hicsparse::type();

//----------------------------------------------------------------------------------------------------------------------

template<typename ValueTy>
atlas::linalg::SparseMatrixStorage createSparseMatrixStorage() {
    auto S = make_sparse_matrix_storage<ValueTy>(
        eckit::linalg::SparseMatrix{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}}
    );
    S.updateDevice();
    return S;
}

//----------------------------------------------------------------------------------------------------------------------

CASE("sparse-matrix vector multiply (spmv) in double [backend=hicsparse]") {
    // "square" matrix
    // A =  2  . -3
    //      .  2  .
    //      .  .  2
    // x = 1 2 3
    // y = 1 2 3
    
    sparse::current_backend(hicsparse);
    
    auto A = createSparseMatrixStorage<double>();
    auto A_device_view = make_device_view<double,int>(A);
    
    SECTION("View of atlas::Array [backend=hicsparse]") {
        ArrayVector<double> x(Vector{1., 2., 3.});
        ArrayVector<double> y(3);
        const auto x_device_view = x.device_view();
        auto y_device_view = y.device_view();
        sparse_matrix_multiply(A_device_view, x_device_view, y_device_view);
        expect_equal(y.host_view(), ArrayVector<double>(Vector{-7., 4., 6.}).host_view());
        // sparse_matrix_multiply of sparse matrix and vector of non-matching sizes should fail
        {
            ArrayVector<double> x2(2);
            auto x2_device_view = x2.device_view();
            EXPECT_THROWS_AS(sparse_matrix_multiply(A_device_view, x2_device_view, y_device_view), eckit::AssertionFailed);
        }
    }

    SECTION("sparse_matrix_multiply_add [backend=hicsparse]") {
        ArrayVector<double> x(Vector{1., 2., 3.});
        ArrayVector<double> y(Vector{4., 5., 6.});
        auto x_device_view = x.device_view();
        auto y_device_view = y.device_view();
        sparse_matrix_multiply_add(A_device_view, x_device_view, y_device_view);
        expect_equal(y.host_view(), ArrayVector<double>(Vector{-3., 9., 12.}).host_view());
        // sparse_matrix_multiply of sparse matrix and vector of non-matching sizes should fail
        {
            ArrayVector<double> x2(2);
            auto x2_device_view = x2.device_view();
            EXPECT_THROWS_AS(sparse_matrix_multiply_add(A_device_view, x2_device_view, y_device_view), eckit::AssertionFailed);
        }
    }
}

CASE("sparse-matrix vector multiply (spmv) in float [backend=hicsparse]") {
    // "square" matrix
    // A =  2  . -3
    //      .  2  .
    //      .  .  2
    // x = 1 2 3
    // y = 1 2 3
    
    sparse::current_backend(hicsparse);
    
    auto A = createSparseMatrixStorage<float>();
    auto A_device_view = make_device_view<float,int>(A);
    
    SECTION("View of atlas::Array [backend=hicsparse]") {
        ArrayVector<float> x(Vector{1., 2., 3.});
        ArrayVector<float> y(3);
        const auto x_device_view = x.device_view();
        auto y_device_view = y.device_view();
        sparse_matrix_multiply(A_device_view, x_device_view, y_device_view);
        expect_equal(y.host_view(), ArrayVector<float>(Vector{-7., 4., 6.}).host_view());
        // sparse_matrix_multiply of sparse matrix and vector of non-matching sizes should fail
        {
            ArrayVector<float> x2(2);
            auto x2_device_view = x2.device_view();
            EXPECT_THROWS_AS(sparse_matrix_multiply(A_device_view, x2_device_view, y_device_view), eckit::AssertionFailed);
        }
    }
}



CASE("sparse-matrix matrix multiply (spmm) in double [backend=hicsparse]") {
    // "square"
    // A =  2  . -3
    //      .  2  .
    //      .  .  2
    // x = 1 2 3
    // y = 1 2 3
    sparse::current_backend(hicsparse);

    auto A = createSparseMatrixStorage<double>();
    auto A_device_view = make_device_view<double,int>(A);
    Matrix m{{1., 2.}, {3., 4.}, {5., 6.}};
    Matrix c_exp{{-13., -14.}, {6., 8.}, {10., 12.}};

    SECTION("View of atlas::Array PointsRight [backend=hicsparse]") {
        ArrayMatrix<double, Indexing::layout_right> ma(m);
        ArrayMatrix<double, Indexing::layout_right> c(3, 2);
        auto ma_device_view = ma.device_view();
        auto c_device_view = c.device_view();
        sparse_matrix_multiply(A_device_view, ma_device_view, c_device_view, Indexing::layout_right);
        expect_equal(c.host_view(), ArrayMatrix<double, Indexing::layout_right>(c_exp).host_view());
    }

    SECTION("sparse_matrix_multiply [backend=hicsparse]") {
        auto backend = sparse::backend::hicsparse();
        ArrayMatrix<double> ma(m);
        ArrayMatrix<double> c(3, 2);
        auto ma_device_view = ma.device_view();
        auto c_device_view = c.device_view();
        sparse_matrix_multiply(A_device_view, ma_device_view, c_device_view, backend);
        expect_equal(c.host_view(), ArrayMatrix<double>(c_exp).host_view());
    }

    SECTION("SparseMatrixMultiply [backend=hicsparse] 1") {
        auto spmm = SparseMatrixMultiply{sparse::backend::hicsparse()};
        ArrayMatrix<double> ma(m);
        ArrayMatrix<double> c(3, 2);
        auto ma_device_view = ma.device_view();
        auto c_device_view = c.device_view();
        spmm(A_device_view, ma_device_view, c_device_view);
        expect_equal(c.host_view(), ArrayMatrix<double>(c_exp).host_view());
    }

    SECTION("SparseMatrixMultiply [backend=hicsparse] 2") {
        auto spmm = SparseMatrixMultiply{hicsparse};
        ArrayMatrix<double> ma(m);
        ArrayMatrix<double> c(3, 2);
        auto ma_device_view = ma.device_view();
        auto c_device_view = c.device_view();
        spmm(A_device_view, ma_device_view, c_device_view);
        expect_equal(c.host_view(), ArrayMatrix<double>(c_exp).host_view());
    }

    SECTION("sparse_matrix_multiply_add [backend=hicsparse]") {
        ArrayMatrix<double> x(m);
        ArrayMatrix<double> y(m);
        Matrix y_exp{{-12., -12.}, {9., 12.}, {15., 18.}};
        auto x_device_view = x.device_view();
        auto y_device_view = y.device_view();
        sparse_matrix_multiply_add(A_device_view, x_device_view, y_device_view, sparse::backend::hicsparse());
        expect_equal(y.host_view(), ArrayMatrix<double>(y_exp).host_view());
    }
}

CASE("sparse-matrix matrix multiply (spmm) in float [backend=hicsparse]") {
    // "square"
    // A =  2  . -3
    //      .  2  .
    //      .  .  2
    // x = 1 2 3
    // y = 1 2 3
    sparse::current_backend(hicsparse);

    auto A = createSparseMatrixStorage<float>();
    auto A_device_view = make_device_view<float,int>(A);
    Matrix m{{1., 2.}, {3., 4.}, {5., 6.}};
    Matrix c_exp{{-13., -14.}, {6., 8.}, {10., 12.}};

    SECTION("View of atlas::Array PointsRight [backend=hicsparse]") {
        ArrayMatrix<float, Indexing::layout_right> ma(m);
        ArrayMatrix<float, Indexing::layout_right> c(3, 2);
        auto ma_device_view = ma.device_view();
        auto c_device_view = c.device_view();
        sparse_matrix_multiply(A_device_view, ma_device_view, c_device_view, Indexing::layout_right);
        expect_equal(c.host_view(), ArrayMatrix<float, Indexing::layout_right>(c_exp).host_view());
    }
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}