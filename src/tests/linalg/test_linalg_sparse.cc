/*
 * (C) Copyright 1996- ECMWF.
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
static std::string eckit_linalg = sparse::backend::eckit_linalg::type();
static std::string openmp       = sparse::backend::openmp::type();

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
    array::ArrayView<Value, 2>& view() { return view_; }
    ArrayMatrix(const eckit::linalg::Matrix& m): ArrayMatrix(m.rows(), m.cols()) {
        for (int r = 0; r < m.rows(); ++r) {
            for (int c = 0; c < m.cols(); ++c) {
                auto& v = layout_left ? view_(r, c) : view_(c, r);
                v       = m(r, c);
            }
        }
    }
    ArrayMatrix(int r, int c): array(make_shape(r, c)), view_(array::make_view<Value, 2>(array)) {}

private:
    static constexpr bool layout_left = (indexing == Indexing::layout_left);
    static array::ArrayShape make_shape(int rows, int cols) {
        return layout_left ? array::make_shape(rows, cols) : array::make_shape(cols, rows);
    }
    array::ArrayT<Value> array;
    array::ArrayView<Value, 2> view_;
};

// 1D array constructable from eckit::linalg::Vector
template <typename Value>
struct ArrayVector {
    array::ArrayView<Value, 1>& view() { return view_; }
    ArrayVector(const eckit::linalg::Vector& v): ArrayVector(v.size()) {
        for (int n = 0; n < v.size(); ++n) {
            view_[n] = v[n];
        }
    }
    ArrayVector(int size): array(size), view_(array::make_view<Value, 1>(array)) {}

private:
    array::ArrayT<Value> array;
    array::ArrayView<Value, 1> view_;
};

bool operator==(const SparseMatrix& A, const SparseMatrix& B) {
    if (A.rows() != B.rows() || A.cols() != B.cols() || A.nonZeros() != B.nonZeros()) {
        return false;
    }
    const auto A_data_view = eckit::testing::make_view(A.data(), A.nonZeros());
    const auto A_outer_view = eckit::testing::make_view(A.outer(), A.rows()+1);
    const auto A_inner_view = eckit::testing::make_view(A.inner(), A.nonZeros());
    const auto B_data_view = eckit::testing::make_view(B.data(), B.nonZeros());
    const auto B_outer_view = eckit::testing::make_view(B.outer(), B.rows()+1);
    const auto B_inner_view = eckit::testing::make_view(B.inner(), B.nonZeros());
    if (A_data_view != B_data_view ||
        A_outer_view != B_outer_view ||
        A_inner_view != B_inner_view) {
        return false;
    }
    return true;
}

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

    sparse::default_backend(eckit_linalg).set("backend", "generic");
    EXPECT_EQ(sparse::default_backend(eckit_linalg).getString("backend"), "generic");

    const sparse::Backend backend_default      = sparse::Backend();
    const sparse::Backend backend_openmp       = sparse::backend::openmp();
    const sparse::Backend backend_eckit_linalg = sparse::backend::eckit_linalg();
    EXPECT_EQ(backend_default.type(), openmp);
    EXPECT_EQ(backend_openmp.type(), openmp);
    EXPECT_EQ(backend_eckit_linalg.type(), eckit_linalg);

    EXPECT_EQ(std::string(backend_openmp), openmp);
    EXPECT_EQ(std::string(backend_eckit_linalg), eckit_linalg);
}

//----------------------------------------------------------------------------------------------------------------------

CASE("SparseMatrix default constructor") {
    SparseMatrix A;
    EXPECT_EQ(A.rows(), 0);
    EXPECT_EQ(A.cols(), 0);
    EXPECT_EQ(A.nonZeros(), 0);
}

CASE("SparseMatrix copy constructor") {
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    SparseMatrix B{A};
    EXPECT(A == B);
}

CASE("SparseMatrix assignment constructor") {
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    auto B = A;
    EXPECT(A == B);
}

CASE("SparseMatrix assignment") {
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    SparseMatrix B;
    B = A;
    EXPECT(A == B);
}

CASE("SparseMatrix triplet constructor") {
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    const auto A_data_view = eckit::testing::make_view(A.data(), A.nonZeros());
    const auto A_outer_view = eckit::testing::make_view(A.outer(), A.rows()+1);
    const auto A_inner_view = eckit::testing::make_view(A.inner(), A.nonZeros());
    
    EXPECT_EQ(A.rows(), 3);
    EXPECT_EQ(A.cols(), 3);
    EXPECT_EQ(A.nonZeros(), 4);
    
    const std::vector<SparseMatrix::Scalar> test_data{2., -3., 2., 2.};
    const std::vector<SparseMatrix::Index> test_outer{0, 2, 3, 4};
    const std::vector<SparseMatrix::Index> test_inner{0, 2, 1, 2};
    const auto test_data_view = eckit::testing::make_view(test_data.data(), test_data.size());
    const auto test_outer_view = eckit::testing::make_view(test_outer.data(), test_outer.size());
    const auto test_inner_view = eckit::testing::make_view(test_inner.data(), test_inner.size());
    
    EXPECT(A_data_view == test_data_view);
    EXPECT(A_outer_view == test_outer_view);
    EXPECT(A_inner_view == test_inner_view);
}

CASE("SparseMatrix triplet constructor 2") {
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, 0.}, {1, 1, 2.}, {2, 2, 2.}}};
    const auto A_data_view = eckit::testing::make_view(A.data(), A.nonZeros());
    const auto A_outer_view = eckit::testing::make_view(A.outer(), A.rows()+1);
    const auto A_inner_view = eckit::testing::make_view(A.inner(), A.nonZeros());
    
    EXPECT_EQ(A.rows(), 3);
    EXPECT_EQ(A.cols(), 3);
    EXPECT_EQ(A.nonZeros(), 3);
    
    const std::vector<SparseMatrix::Scalar> test_data{2., 2., 2.};
    const std::vector<SparseMatrix::Index> test_outer{0, 1, 2, 3};
    const std::vector<SparseMatrix::Index> test_inner{0, 1, 2};
    const auto test_data_view = eckit::testing::make_view(test_data.data(), test_data.size());
    const auto test_outer_view = eckit::testing::make_view(test_outer.data(), test_outer.size());
    const auto test_inner_view = eckit::testing::make_view(test_inner.data(), test_inner.size());
    
    EXPECT(A_data_view == test_data_view);
    EXPECT(A_outer_view == test_outer_view);
    EXPECT(A_inner_view == test_inner_view);
}

CASE("SparseMatrix swap") {
    SparseMatrix A_test{3, 3, {{0, 0, 2.}, {0, 2, 0.}, {1, 1, 2.}, {2, 2, 2.}}};
    SparseMatrix A{A_test};
    SparseMatrix B_test{1, 1, {{0, 0, 7.}}};
    SparseMatrix B{B_test};
    A.swap(B);
    EXPECT(A == B_test);
    EXPECT(B == A_test);
}

CASE("SparseMatrix transpose") {
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    SparseMatrix AT{3, 3, {{0, 0, 2.}, {1, 1, 2.}, {2, 0, -3.}, {2, 2, 2.}}};
    A.transpose();
    EXPECT(A == AT);
}

CASE("SparseMatrix prune") {
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, 0}, {1, 1, 2.}, {2, 2, 2.}}};
    SparseMatrix A_pruned{3, 3, {{0, 0, 2.}, {1, 1, 2.}, {2, 2, 2.}}};
    A.prune();
    EXPECT(A == A_pruned);
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
            sparse_matrix_multiply(A, x.view(), y.view());
            expect_equal(y.view(), Vector{-7., 4., 6.});
            // sparse_matrix_multiply of sparse matrix and vector of non-matching sizes should fail
            {
                ArrayVector<double> x2(2);
                EXPECT_THROWS_AS(sparse_matrix_multiply(A, x2.view(), y.view()), eckit::AssertionFailed);
            }
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
            sparse_matrix_multiply(A, ma.view(), c.view(), Indexing::layout_right);
            expect_equal(c.view(), c_exp);
        }
    }

    SECTION("sparse_matrix_multiply [backend=openmp]") {
        sparse::current_backend(eckit_linalg);  // expected to be ignored further
        auto backend = sparse::backend::openmp();
        ArrayMatrix<float> ma(m);
        ArrayMatrix<float> c(3, 2);
        sparse_matrix_multiply(A, ma.view(), c.view(), backend);
        expect_equal(c.view(), ArrayMatrix<float>(c_exp).view());
    }

    SECTION("SparseMatrixMultiply [backend=openmp] 1") {
        sparse::current_backend(eckit_linalg);  // expected to be ignored
        auto spmm = SparseMatrixMultiply{sparse::backend::openmp()};
        ArrayMatrix<float> ma(m);
        ArrayMatrix<float> c(3, 2);
        spmm(A, ma.view(), c.view());
        expect_equal(c.view(), ArrayMatrix<float>(c_exp).view());
    }

    SECTION("SparseMatrixMultiply [backend=openmp] 2") {
        sparse::current_backend(eckit_linalg);  // expected to be ignored
        auto spmm = SparseMatrixMultiply{openmp};
        ArrayMatrix<float> ma(m);
        ArrayMatrix<float> c(3, 2);
        spmm(A, ma.view(), c.view());
        expect_equal(c.view(), ArrayMatrix<float>(c_exp).view());
    }
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
