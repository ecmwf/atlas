/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

// #include "atlas/linalg/SparseMatrix.h"

#include "hic/hic.h"
#include "hic/hicsparse.h"

#include "tests/AtlasTestEnvironment.h"


#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/linalg/sparse/SparseMatrixView.h"
#include "atlas/linalg/sparse/MakeSparseMatrixStorageEigen.h"
#include "atlas/linalg/sparse/MakeSparseMatrixStorageEckit.h"

using atlas::linalg::SparseMatrixStorage;
using atlas::linalg::make_sparse_matrix_storage;
using atlas::linalg::make_host_view;
using atlas::linalg::make_device_view;

namespace atlas {
namespace test {

template<typename MatrixValue, typename IndexType, typename XValue, typename YValue>
void do_hicsparse_matrix_multiply(const atlas::linalg::SparseMatrixView<MatrixValue,IndexType>& device_A, const XValue* device_x, YValue* device_y) {
    // Create sparse library handle
    hicsparseHandle_t handle;
    HICSPARSE_CALL( hicsparseCreate(&handle) );

    auto get_hic_value_type = [](const auto& dummy) -> hicDataType {
        if(std::is_same_v<std::decay_t<decltype(dummy)>,double> ) {
            return HIC_R_64F;
        }
        else if (std::is_same_v<std::decay_t<decltype(dummy)>,float> ) {
            return HIC_R_32F;
        }
        else if (std::is_same_v<std::decay_t<decltype(dummy)>,int> ) {
            return HIC_R_32F;
        }
        ATLAS_NOTIMPLEMENTED;
    };
    auto get_hic_index_type = [](const auto& dummy) -> hicsparseIndexType_t {
        if (std::is_same_v<std::make_signed_t<std::decay_t<decltype(dummy)>>,int> ) {
            return HICSPARSE_INDEX_32I;
        }
        else if (std::is_same_v<std::make_signed_t<std::decay_t<decltype(dummy)>>,long> ) {
            return HICSPARSE_INDEX_64I;
        }
        ATLAS_NOTIMPLEMENTED;
    };
    auto compute_type = get_hic_value_type(MatrixValue());
    auto index_type   = get_hic_index_type(IndexType());

    hicsparseConstSpMatDescr_t matA;
    HICSPARSE_CALL( hicsparseCreateConstCsr(
        &matA,
        device_A.rows(), device_A.cols(), device_A.nnz(),
        device_A.outer(),    // row_offsets
        device_A.inner(),    // column_indices
        device_A.value(),    // values
        index_type,
        index_type,
        HICSPARSE_INDEX_BASE_ZERO,
        get_hic_value_type(MatrixValue())) );

    // Create dense matrix descriptors
    hicsparseConstDnVecDescr_t vecX;
    HICSPARSE_CALL( hicsparseCreateConstDnVec(
        &vecX,
        device_A.cols(),
        device_x,
        get_hic_value_type(XValue())) );

    hicsparseDnVecDescr_t vecY;
    HICSPARSE_CALL( hicsparseCreateDnVec(
        &vecY,
        device_A.rows(),
        device_y,
        get_hic_value_type(YValue())) );

    // Set parameters in compute_type
    const MatrixValue alpha = 1.0;
    const MatrixValue beta = 0.0;

    // Determine buffer size
    size_t bufferSize = 0;
    HICSPARSE_CALL( hicsparseSpMV_bufferSize(
        handle,
        HICSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha,
        matA,
        vecX,
        &beta,
        vecY,
        compute_type,
        HICSPARSE_SPMV_ALG_DEFAULT,
        &bufferSize) );

    // Allocate buffer
    char* buffer;
    HIC_CALL( hicMalloc(&buffer, bufferSize) );

    // Perform SpMV
    HICSPARSE_CALL( hicsparseSpMV(
        handle,
        HICSPARSE_OPERATION_NON_TRANSPOSE,
        &alpha,
        matA,
        vecX,
        &beta,
        vecY,
        compute_type,
        HICSPARSE_SPMV_ALG_DEFAULT,
        buffer) );

    HIC_CALL( hicFree(buffer) );
    HICSPARSE_CALL( hicsparseDestroyDnVec(vecY) );
    HICSPARSE_CALL( hicsparseDestroyDnVec(vecX) );
    HICSPARSE_CALL( hicsparseDestroySpMat(matA) );
    HICSPARSE_CALL( hicsparseDestroy(handle) );
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

    m.updateDevice();
    auto host_matrix_view   = make_host_view<Value,Index>(m);
    auto device_matrix_view = make_device_view<Value,Index>(m);

    std::vector<Index>  host_outer(host_matrix_view.outer_size());
    std::vector<Index>  host_inner(host_matrix_view.inner_size());
    std::vector<Value>  host_value(host_matrix_view.value_size());

    hicMemcpy(host_outer.data(), device_matrix_view.outer(), device_matrix_view.outer_size() * sizeof(Index),  hicMemcpyDeviceToHost);
    hicMemcpy(host_inner.data(), device_matrix_view.inner(), device_matrix_view.inner_size() * sizeof(Index),  hicMemcpyDeviceToHost);
    hicMemcpy(host_value.data(), device_matrix_view.value(), device_matrix_view.value_size() * sizeof(Value),  hicMemcpyDeviceToHost);

    EXPECT(host_outer == expected_outer);
    EXPECT(host_inner == expected_inner);
    EXPECT(host_value == expected_value);

    array::ArrayT<Value> x(m.cols());
    array::ArrayT<Value> y(m.rows());
    array::make_host_view<Value,1>(x).assign({1.0, 2.0, 3.0});

    m.updateDevice();
    x.updateDevice();
    y.allocateDevice();

    do_hicsparse_matrix_multiply(device_matrix_view, x.device_data(), y.device_data());

    // Check result
    y.updateHost();
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
    EXPECT_THROWS(make_host_view<double>(S));        // value data type mismatch
    EXPECT_THROWS((make_host_view<float,long>(S)));  // index data type mismatch
}
CASE("Test exceptions in make_device_view") {
    SparseMatrixStorage S( create_sparsematrix_storage<float>() );
    EXPECT_THROWS(make_device_view<float>(S)); // device is not allocated
    S.updateDevice();
    EXPECT_THROWS(make_device_view<double>(S));       // value data type mismatch
    EXPECT_THROWS((make_device_view<float,long>(S))); // index data type mismatch
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
