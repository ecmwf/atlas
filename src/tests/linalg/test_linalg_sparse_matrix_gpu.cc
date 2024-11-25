/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/linalg/SparseMatrix.h"

#include "hic/hic.h"
#include "hic/hicsparse.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::linalg;

namespace atlas {
namespace test {

CASE("SparseMatrix update device") {
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    
    A.updateDevice();
    
    std::vector<SparseMatrix::Scalar> A_data_cpy(A.data_size());
    std::vector<SparseMatrix::Index>  A_outer_cpy(A.outer_size());
    std::vector<SparseMatrix::Index>  A_inner_cpy(A.inner_size());

    hicMemcpy(A_data_cpy.data(),  A.device_data(),  A.data_size()  * sizeof(SparseMatrix::Scalar), hicMemcpyDeviceToHost);
    hicMemcpy(A_outer_cpy.data(), A.device_outer(), A.outer_size() * sizeof(SparseMatrix::Index),  hicMemcpyDeviceToHost);
    hicMemcpy(A_inner_cpy.data(), A.device_inner(), A.inner_size() * sizeof(SparseMatrix::Index),  hicMemcpyDeviceToHost);

    std::vector<SparseMatrix::Scalar> A_data_test{2., -3., 2., 2.};
    std::vector<SparseMatrix::Index>  A_outer_test{0, 2, 3, 4};
    std::vector<SparseMatrix::Index>  A_inner_test{0, 2, 1, 2};

    EXPECT(A_data_cpy == A_data_test);
    EXPECT(A_outer_cpy == A_outer_test);
    EXPECT(A_inner_cpy == A_inner_test);
}

void do_hicsparse_matrix_multiply(const SparseMatrix& A, const int N, const double* device_x, double* device_y);

CASE("SparseMatrix hicsparse multiply") {
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    array::ArrayT<double> x(3);
    array::ArrayT<double> y(3);
    array::make_host_view<double,1>(x).assign({1.0, 2.0, 3.0});

    A.updateDevice();
    x.updateDevice();
    y.allocateDevice();

    do_hicsparse_matrix_multiply(
            A,
            x.shape(0),
            x.device_data<double>(),
            y.device_data<double>());

    // Check result
    constexpr int N=3;
    const double expected_y[N] = {-7.0, 4.0, 6.0};
    y.updateHost();
    auto hy = array::make_host_view<double,1>(y);
    for( int i = 0; i < N; ++i ) {
        EXPECT_EQ(hy[i], expected_y[i]);
    }
}

void do_hicsparse_matrix_multiply(const SparseMatrix& A, const int N, const double* device_x, double* device_y) {
    // Create sparse library handle
    hicsparseHandle_t handle;
    HICSPARSE_CALL( hicsparseCreate(&handle) );

    static_assert(std::is_same_v<SparseMatrix::Index, int>,"Following assumes Index==int");
    static_assert(std::is_same_v<SparseMatrix::Scalar, double>,"Following assumes Index==int");
    hicsparseConstSpMatDescr_t matA;
    HICSPARSE_CALL( hicsparseCreateConstCsr(
        &matA,
        A.rows(), A.cols(), A.nonZeros(),
        A.device_outer(),    // row_offsets
        A.device_inner(),    // column_indices
        A.device_data(),     // values
        HICSPARSE_INDEX_32I,
        HICSPARSE_INDEX_32I,
        HICSPARSE_INDEX_BASE_ZERO,
        HIC_R_64F) );

    // Create dense matrix descriptors
    hicsparseConstDnVecDescr_t vecX;
    HICSPARSE_CALL( hicsparseCreateConstDnVec(
        &vecX,
        N,
        device_x,
        HIC_R_64F) );

    hicsparseDnVecDescr_t vecY;
    HICSPARSE_CALL( hicsparseCreateDnVec(
        &vecY,
        N,
        device_y,
        HIC_R_64F) );

    // Set parameters
    const double alpha = 1.0;
    const double beta = 0.0;

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
        HIC_R_64F,
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
        HIC_R_64F,
        HICSPARSE_SPMV_ALG_DEFAULT,
        buffer) );

    HIC_CALL( hicFree(buffer) );
    HICSPARSE_CALL( hicsparseDestroyDnVec(vecY) );
    HICSPARSE_CALL( hicsparseDestroyDnVec(vecX) );
    HICSPARSE_CALL( hicsparseDestroySpMat(matA) );
    HICSPARSE_CALL( hicsparseDestroy(handle) );
}


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
