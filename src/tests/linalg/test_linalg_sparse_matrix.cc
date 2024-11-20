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

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::linalg;

namespace atlas {
namespace test {

bool operator==(const SparseMatrix& A, const SparseMatrix& B) {
    if (A.rows() != B.rows() || A.cols() != B.cols() || A.nonZeros() != B.nonZeros()) {
        return false;
    }
    const auto A_data_view = eckit::testing::make_view(A.data(), A.data_size());
    const auto A_outer_view = eckit::testing::make_view(A.outer(), A.outer_size());
    const auto A_inner_view = eckit::testing::make_view(A.inner(), A.inner_size());
    const auto B_data_view = eckit::testing::make_view(B.data(), B.data_size());
    const auto B_outer_view = eckit::testing::make_view(B.outer(), B.outer_size());
    const auto B_inner_view = eckit::testing::make_view(B.inner(), B.inner_size());
    if (A_data_view != B_data_view ||
        A_outer_view != B_outer_view ||
        A_inner_view != B_inner_view) {
        return false;
    }
    return true;
}

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

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
