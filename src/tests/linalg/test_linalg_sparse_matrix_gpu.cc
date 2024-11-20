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

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::linalg;

namespace atlas {
namespace test {

CASE("SparseMatrix update device") {
    SparseMatrix A{3, 3, {{0, 0, 2.}, {0, 2, -3.}, {1, 1, 2.}, {2, 2, 2.}}};
    
    A.updateDevice();
    
    std::vector<SparseMatrix::Scalar> A_data_cpy(A.data_size());
    std::vector<SparseMatrix::Index> A_outer_cpy(A.outer_size());
    std::vector<SparseMatrix::Index> A_inner_cpy(A.inner_size());

    hicMemcpy(A_data_cpy.data(), A.device_data(), A.data_size() * sizeof(SparseMatrix::Scalar), hicMemcpyDeviceToHost);
    hicMemcpy(A_outer_cpy.data(), A.device_outer(), A.outer_size() * sizeof(SparseMatrix::Index), hicMemcpyDeviceToHost);
    hicMemcpy(A_inner_cpy.data(), A.device_inner(), A.inner_size() * sizeof(SparseMatrix::Index), hicMemcpyDeviceToHost);

    std::vector<SparseMatrix::Scalar> A_data_test{2., -3., 2., 2.};
    std::vector<SparseMatrix::Index> A_outer_test{0, 2, 3, 4};
    std::vector<SparseMatrix::Index> A_inner_test{0, 2, 1, 2};

    EXPECT(A_data_cpy == A_data_test);
    EXPECT(A_outer_cpy == A_outer_test);
    EXPECT(A_inner_cpy == A_inner_test);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
