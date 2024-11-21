/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <tuple>
#include <vector>

#include "eckit/linalg/Matrix.h"
#include "eckit/linalg/Vector.h"

#include "atlas/array.h"
#include "atlas/linalg/sparse.h"

#include "tests/AtlasTestEnvironment.h"


using namespace atlas::linalg;

namespace atlas {
namespace test {

//----------------------------------------------------------------------------------------------------------------------

// strings to be used in the tests
static std::string hicsparse    = sparse::backend::hicsparse::type();

//----------------------------------------------------------------------------------------------------------------------

// Only reason to define these derived classes is for nicer constructors and convenience in the tests

class Vector : public eckit::linalg::Vector {
public:
    using Scalar = eckit::linalg::Scalar;
    using eckit::linalg::Vector::Vector;
    Vector(const std::initializer_list<Scalar>& v): eckit::linalg::Vector::Vector(v.size()) {
        size_t i = 0;
        for (auto& s : v) {
            operator[](i++) = s;
        }
    }
};

class Matrix : public eckit::linalg::Matrix {
public:
    using Scalar = eckit::linalg::Scalar;
    using eckit::linalg::Matrix::Matrix;

    Matrix(const std::initializer_list<std::vector<Scalar>>& m):
        eckit::linalg::Matrix::Matrix(m.size(), m.size() ? m.begin()->size() : 0) {
        size_t r = 0;
        for (auto& row : m) {
            for (size_t c = 0; c < cols(); ++c) {
                operator()(r, c) = row[c];
            }
            ++r;
        }
    }
};

// 2D array constructable from eckit::linalg::Matrix
// Indexing/memorylayout and data type can be customized for testing
template <typename Value, Indexing indexing = Indexing::layout_left>
struct ArrayMatrix {
    array::ArrayView<Value, 2> view() {
        array.syncHostDevice();
        return array::make_view<Value, 2>(array);
    }
    array::ArrayView<Value, 2> device_view() {
        array.syncHostDevice();
        return array::make_device_view<Value, 2>(array);
    }
    void setHostNeedsUpdate(bool b) {
        array.setHostNeedsUpdate(b);
    }
    ArrayMatrix(const eckit::linalg::Matrix& m): ArrayMatrix(m.rows(), m.cols()) {
        auto view_ = array::make_view<Value, 2>(array);
        for (int r = 0; r < m.rows(); ++r) {
            for (int c = 0; c < m.cols(); ++c) {
                auto& v = layout_left ? view_(r, c) : view_(c, r);
                v       = m(r, c);
            }
        }
    }
    ArrayMatrix(int r, int c): array(make_shape(r, c)) {}

private:
    static constexpr bool layout_left = (indexing == Indexing::layout_left);
    static array::ArrayShape make_shape(int rows, int cols) {
        return layout_left ? array::make_shape(rows, cols) : array::make_shape(cols, rows);
    }
    array::ArrayT<Value> array;
};

// 1D array constructable from eckit::linalg::Vector
template <typename Value>
struct ArrayVector {
    array::ArrayView<Value, 1> view() {
        array.syncHostDevice();
        return array::make_view<Value, 1>(array);
    }
    array::ArrayView<const Value, 1> const_view() {
        array.syncHostDevice();
        return array::make_view<const Value, 1>(array);
    }
    array::ArrayView<Value, 1> device_view() {
        array.syncHostDevice();
        return array::make_device_view<Value, 1>(array);
    }
    void setHostNeedsUpdate(bool b) {
        array.setHostNeedsUpdate(b);
    }
    ArrayVector(const eckit::linalg::Vector& v): ArrayVector(v.size()) {
        auto view_ = array::make_view<Value, 1>(array);
        for (int n = 0; n < v.size(); ++n) {
            view_[n] = v[n];
        }
    }
    ArrayVector(int size): array(size) {}

private:
    array::ArrayT<Value> array;
};

//----------------------------------------------------------------------------------------------------------------------

template <typename T>
void expect_equal(T* v, T* r, size_t s) {
    EXPECT(is_approximately_equal(eckit::testing::make_view(v, s), eckit::testing::make_view(r, s), T(1.e-5)));
}

template <class T1, class T2>
void expect_equal(const T1& v, const T2& r) {
    expect_equal(v.data(), r.data(), std::min(v.size(), r.size()));
}

//----------------------------------------------------------------------------------------------------------------------

CASE("sparse-matrix vector multiply (spmv) [backend=hicsparse]") {
    // "square" matrix
    // A =  2  . -3
    //      .  2  .
    //      .  .  2
    // x = 1 2 3
    // y = 1 2 3
    
    sparse::current_backend(hicsparse);
    
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    
    SECTION("View of atlas::Array [backend=hicsparse]") {
        ArrayVector<double> x(Vector{1., 2., 3.});
        ArrayVector<double> y(3);
        const auto x_device_view = x.device_view();
        auto y_device_view = y.device_view();
        sparse_matrix_multiply(A, x_device_view, y_device_view);
        y.setHostNeedsUpdate(true);
        auto y_view = y.view();
        expect_equal(y.view(), Vector{-7., 4., 6.});
        // sparse_matrix_multiply of sparse matrix and vector of non-matching sizes should fail
        {
            ArrayVector<double> x2(2);
            auto x2_device_view = x2.device_view();
            EXPECT_THROWS_AS(sparse_matrix_multiply(A, x2_device_view, y_device_view), eckit::AssertionFailed);
        }
    }

    SECTION("sparse_matrix_multiply_add [backend=hicsparse]") {
        ArrayVector<double> x(Vector{1., 2., 3.});
        ArrayVector<double> y(Vector{4., 5., 6.});
        auto x_device_view = x.device_view();
        auto y_device_view = y.device_view();
        sparse_matrix_multiply_add(A, x_device_view, y_device_view);
        y.setHostNeedsUpdate(true);
        auto y_view = y.view();
        expect_equal(y.view(), Vector{-3., 9., 12.});
        // sparse_matrix_multiply of sparse matrix and vector of non-matching sizes should fail
        {
            ArrayVector<double> x2(2);
            auto x2_device_view = x2.device_view();
            EXPECT_THROWS_AS(sparse_matrix_multiply_add(A, x2_device_view, y_device_view), eckit::AssertionFailed);
        }
    }
}

CASE("sparse-matrix matrix multiply (spmm) [backend=hicsparse]") {
    // "square"
    // A =  2  . -3
    //      .  2  .
    //      .  .  2
    // x = 1 2 3
    // y = 1 2 3
    sparse::current_backend(hicsparse);

    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    Matrix m{{1., 2.}, {3., 4.}, {5., 6.}};
    Matrix c_exp{{-13., -14.}, {6., 8.}, {10., 12.}};

    SECTION("View of atlas::Array PointsRight [backend=hicsparse]") {
        ArrayMatrix<double, Indexing::layout_right> ma(m);
        ArrayMatrix<double, Indexing::layout_right> c(3, 2);
        auto ma_device_view = ma.device_view();
        auto c_device_view = c.device_view();
        sparse_matrix_multiply(A, ma_device_view, c_device_view, Indexing::layout_right);
        c.setHostNeedsUpdate(true);
        auto c_view = c.view();
        expect_equal(c_view, ArrayMatrix<double, Indexing::layout_right>(c_exp).view());
    }

    SECTION("sparse_matrix_multiply [backend=hicsparse]") {
        auto backend = sparse::backend::hicsparse();
        ArrayMatrix<double> ma(m);
        ArrayMatrix<double> c(3, 2);
        auto ma_device_view = ma.device_view();
        auto c_device_view = c.device_view();
        sparse_matrix_multiply(A, ma_device_view, c_device_view, backend);
        c.setHostNeedsUpdate(true);
        auto c_view = c.view();
        expect_equal(c_view, ArrayMatrix<double>(c_exp).view());
    }

    SECTION("SparseMatrixMultiply [backend=hicsparse] 1") {
        auto spmm = SparseMatrixMultiply{sparse::backend::hicsparse()};
        ArrayMatrix<double> ma(m);
        ArrayMatrix<double> c(3, 2);
        auto ma_device_view = ma.device_view();
        auto c_device_view = c.device_view();
        spmm(A, ma_device_view, c_device_view);
        c.setHostNeedsUpdate(true);
        auto c_view = c.view();
        expect_equal(c_view, ArrayMatrix<double>(c_exp).view());
    }

    SECTION("SparseMatrixMultiply [backend=hicsparse] 2") {
        auto spmm = SparseMatrixMultiply{hicsparse};
        ArrayMatrix<double> ma(m);
        ArrayMatrix<double> c(3, 2);
        auto ma_device_view = ma.device_view();
        auto c_device_view = c.device_view();
        spmm(A, ma_device_view, c_device_view);
        c.setHostNeedsUpdate(true);
        auto c_view = c.view();
        expect_equal(c_view, ArrayMatrix<double>(c_exp).view());
    }

    SECTION("sparse_matrix_multiply_add [backend=hicsparse]") {
        ArrayMatrix<double> x(m);
        ArrayMatrix<double> y(m);
        Matrix y_exp{{-12., -12.}, {9., 12.}, {15., 18.}};
        auto x_device_view = x.device_view();
        auto y_device_view = y.device_view();
        sparse_matrix_multiply_add(A, x_device_view, y_device_view, sparse::backend::hicsparse());
        y.setHostNeedsUpdate(true);
        auto y_view = y.view();
        expect_equal(y_view, ArrayMatrix<double>(y_exp).view());
    }
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
