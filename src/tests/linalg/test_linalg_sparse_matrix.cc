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


#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/linalg/sparse/SparseMatrixView.h"
#include "atlas/linalg/sparse/MakeSparseMatrixStorageEigen.h"
#include "atlas/linalg/sparse/MakeSparseMatrixStorageEckit.h"

using atlas::linalg::SparseMatrixStorage;
using atlas::linalg::make_sparse_matrix_storage;

namespace atlas {
namespace test {

template<typename MatrixValue, typename IndexType, typename XValue, typename YValue>
void do_matrix_multiply(const atlas::linalg::SparseMatrixView<MatrixValue,IndexType>& W, const XValue* x, YValue* y) {  
    const auto outer  = W.outer();
    const auto inner  = W.inner();
    const auto weight = W.value();
    int rows          = W.rows();
    for (int r = 0; r < rows; ++r) {
        y[r] = 0.;
        for (IndexType c = outer[r]; c < outer[r + 1]; ++c) {
            IndexType n = inner[c];
            y[r] += weight[c] * x[n];
        }
    }
}

#if ATLAS_HAVE_EIGEN
template <typename Value, typename Index = std::make_signed_t<eckit::linalg::Index>>
Eigen::SparseMatrix<Value, Eigen::RowMajor, Index> create_eigen_sparsematrix() {
    Eigen::SparseMatrix<Value, Eigen::RowMajor, Index> matrix(3,3);
    std::vector<Eigen::Triplet<Value>> triplets;
    triplets.reserve(4);
    triplets.emplace_back(0, 0,  2.);
    triplets.emplace_back(0, 2, -3.);
    triplets.emplace_back(1, 1,  2.);
    triplets.emplace_back(2, 2,  2.);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}
#endif

eckit::linalg::SparseMatrix create_eckit_sparsematrix() {
    return eckit::linalg::SparseMatrix{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
}

template <typename Value, typename Index = eckit::linalg::Index>
SparseMatrixStorage create_sparsematrix_storage() {
    return make_sparse_matrix_storage<Value,Index>(create_eckit_sparsematrix());
}

template<typename Value>
void check_array(const array::Array& array, const std::vector<Value>& expected) {
    EXPECT_EQ(array.size(), expected.size());
    auto view = array::make_view<Value,1>(array);
    for( std::size_t i=0; i<expected.size(); ++i) {
        EXPECT_EQ(view[i], expected[i]);
    }
}


template<typename Value, typename Index>
void check_matrix(const SparseMatrixStorage& m) {
    std::vector<Value>  expected_value{2., -3., 2., 2.};
    std::vector<Index>  expected_outer{0, 2, 3, 4};
    std::vector<Index>  expected_inner{0, 2, 1, 2};

    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 3);
    EXPECT_EQ(m.nnz(),  4);

    check_array(m.value(),expected_value);
    check_array(m.outer(),expected_outer);
    check_array(m.inner(),expected_inner);

    auto host_matrix_view   = linalg::make_host_view<Value,Index>(m);

    std::vector<Index>  host_outer(host_matrix_view.outer_size());
    std::vector<Index>  host_inner(host_matrix_view.inner_size());
    std::vector<Value>  host_value(host_matrix_view.value_size());

    array::ArrayT<Value> x(m.cols());
    array::ArrayT<Value> y(m.rows());
    array::make_host_view<Value,1>(x).assign({1.0, 2.0, 3.0});

    do_matrix_multiply(host_matrix_view, x.host_data(), y.host_data());

    // Check result
    std::vector<Value> expected_y{-7.0, 4.0, 6.0};
    check_array(y, expected_y);
}

void check_matrix(const SparseMatrixStorage& m) {
    if (m.value().datatype() == atlas::make_datatype<double>() && m.outer().datatype() == atlas::make_datatype<int>()) {
        check_matrix<double, int>(m);
    }
    else if (m.value().datatype() == atlas::make_datatype<double>() && m.outer().datatype() == atlas::make_datatype<unsigned int>()) {
        check_matrix<double, unsigned int>(m);
    }
    else if (m.value().datatype() == atlas::make_datatype<double>() && m.outer().datatype() == atlas::make_datatype<long>()) {
        check_matrix<double, long>(m);
    }
    else if (m.value().datatype() == atlas::make_datatype<double>() && m.outer().datatype() == atlas::make_datatype<unsigned long>()) {
        check_matrix<double, unsigned long>(m);
    }
    else if (m.value().datatype() == atlas::make_datatype<float>() && m.outer().datatype() == atlas::make_datatype<int>()) {
        check_matrix<float, int>(m);
    }
    else if (m.value().datatype() == atlas::make_datatype<float>() && m.outer().datatype() == atlas::make_datatype<unsigned int>()) {
        check_matrix<float, unsigned int>(m);
    }
    else if (m.value().datatype() == atlas::make_datatype<float>() && m.outer().datatype() == atlas::make_datatype<long>()) {
        check_matrix<float, long>(m);
    }
    else if (m.value().datatype() == atlas::make_datatype<float>() && m.outer().datatype() == atlas::make_datatype<unsigned long>()) {
        check_matrix<float, unsigned long>(m);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

CASE("Create SparseMatrix moving eckit::linalg::SparseMatrix") {
    SparseMatrixStorage S = make_sparse_matrix_storage(create_eckit_sparsematrix());
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Create SparseMatrix moving eckit::linalg::SparseMatrix using std::move") {
    auto A = create_eckit_sparsematrix();
    SparseMatrixStorage S = make_sparse_matrix_storage(std::move(A));
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Create SparseMatrix copying eckit::linalg::SparseMatrix") {
    auto A = create_eckit_sparsematrix();
    SparseMatrixStorage S = make_sparse_matrix_storage(A);
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Create single precision SparseMatrix copying from double precision SparseMatrix ") {
    auto Sdouble = create_sparsematrix_storage<double>();
    SparseMatrixStorage S = make_sparse_matrix_storage<float>(Sdouble);
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Create single precision SparseMatrix moving from double precision SparseMatrix which causes copy") {
    SparseMatrixStorage S = make_sparse_matrix_storage<float>(create_sparsematrix_storage<double>());
    EXPECT_NO_THROW(check_matrix(S));
}
CASE("Create double precision SparseMatrix copying from single precision SparseMatrix ") {
    auto Sfloat = create_sparsematrix_storage<float>();
    SparseMatrixStorage S = make_sparse_matrix_storage<double>(Sfloat);
    EXPECT_NO_THROW(check_matrix(S));
}
CASE("Create double precision SparseMatrix moving from single precision SparseMatrix which causes copy") {
    SparseMatrixStorage S = make_sparse_matrix_storage<double>(create_sparsematrix_storage<float>());
    EXPECT_NO_THROW(check_matrix(S));
}
CASE("Create base with copy constructor from double precision") {
    auto Sdouble = create_sparsematrix_storage<double>();
    SparseMatrixStorage S(Sdouble);
    EXPECT_NO_THROW(check_matrix(S));
}
CASE("Create base with move constructor from double precision") {
    SparseMatrixStorage S(create_sparsematrix_storage<double>());
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Create base with copy constructor from single precision") {
    auto Sfloat = create_sparsematrix_storage<float>();
    SparseMatrixStorage S(Sfloat);
    EXPECT_NO_THROW(check_matrix(S));
}
CASE("Create base with move constructor from single precision") {
    SparseMatrixStorage S(create_sparsematrix_storage<float>());
    EXPECT_NO_THROW(check_matrix(S));
}

#if ATLAS_HAVE_EIGEN
CASE("Copy from Eigen double") {
    auto eigen_matrix = create_eigen_sparsematrix<double>();
    SparseMatrixStorage S = make_sparse_matrix_storage(eigen_matrix);
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Move from Eigen double, avoiding copy") {
    SparseMatrixStorage S = make_sparse_matrix_storage(create_eigen_sparsematrix<double>());
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Copy and convert double from Eigen to single precision") {
    auto eigen_matrix = create_eigen_sparsematrix<double>();
    SparseMatrixStorage S = make_sparse_matrix_storage<float>(eigen_matrix);
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Move and convert double from Eigen to single precision, should cause copy") {
    SparseMatrixStorage S = make_sparse_matrix_storage<float>(create_eigen_sparsematrix<double>());
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Copy from Eigen single") {
    auto eigen_matrix = create_eigen_sparsematrix<float>();
    SparseMatrixStorage S = make_sparse_matrix_storage(eigen_matrix);
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Move from Eigen single, avoiding copy") {
    SparseMatrixStorage S = make_sparse_matrix_storage(create_eigen_sparsematrix<float>());
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Copy and convert single from Eigen to double precision") {
    auto eigen_matrix = create_eigen_sparsematrix<float>();
    SparseMatrixStorage S = make_sparse_matrix_storage<double>(eigen_matrix);
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Move and convert single from Eigen to double precision, should cause copy") {
    SparseMatrixStorage S = make_sparse_matrix_storage<double>(create_eigen_sparsematrix<float>());
    EXPECT_NO_THROW(check_matrix(S));
}

CASE("Create chain of moves from eigen double") {
    SparseMatrixStorage S = make_sparse_matrix_storage<double>(create_eigen_sparsematrix<double>());
    EXPECT_NO_THROW(check_matrix(S));
}
CASE("Create chain of moves from eigen single") {
    SparseMatrixStorage S = make_sparse_matrix_storage<float>(create_eigen_sparsematrix<float>());
    EXPECT_NO_THROW(check_matrix(S));
}
#endif

CASE("Test exceptions in make_host_view") {
    SparseMatrixStorage S( create_sparsematrix_storage<float>() );
    EXPECT_THROWS(linalg::make_host_view<double>(S));        // value data type mismatch
    EXPECT_THROWS((linalg::make_host_view<float,long>(S)));  // index data type mismatch
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
