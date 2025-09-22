/*
 * (C) Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <list>
#include <random>
#include <set>
#include <type_traits>
#include <vector>

#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/linalg/sparse/SparseMatrixTriplet.h"
#include "atlas/runtime/Exception.h"
#include "eckit/testing/Test.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas::test {

using namespace atlas::linalg;

using TripletType = Triplet<double, int>;

CASE("Test triplet.") {
    SECTION("Test triplet construction.") {
        const auto t = Triplet{1, 2, 3.};
        static_assert(std::is_same_v<decltype(t), const TripletType>);
        ATLAS_ASSERT(t.row() == 1);
        ATLAS_ASSERT(t.col() == 2);
        ATLAS_ASSERT(t.value() == 3.);
    }

    SECTION("Test triplet less than operator.") {
        const auto t1 = Triplet{1, 2, 0.};
        const auto t2 = Triplet{1, 3, 0.};
        const auto t3 = Triplet{2, 1, 0.};

        ATLAS_ASSERT(t1 < t2);
        ATLAS_ASSERT(t1 < t3);
        ATLAS_ASSERT(!(t2 < t1));
        ATLAS_ASSERT(!(t3 < t1));
        ATLAS_ASSERT(!(t1 < t1));
    }
}

std::tuple<std::vector<TripletType>, int, int> make_test_triplets() {
    // Return triplet vector, n_rows and n_cols.
    return {{
                {0, 1, 2.},
                {1, 0, 1.},
                {1, 1, 3.},
                {3, 0, 4.},
                {3, 2, 5.},
            },
            4,
            3};
}

CASE("Test sparse matrix construction from triplets.") {
    const auto test_triplets = [](const SparseMatrixStorage& matrix) {
        ATLAS_ASSERT(matrix.rows() == 4);
        ATLAS_ASSERT(matrix.cols() == 3);
        ATLAS_ASSERT(matrix.nnz() == 5);

        const auto outer_view  = array::make_view<int, 1>(matrix.outer());
        const auto inner_view  = array::make_view<int, 1>(matrix.inner());
        const auto values_view = array::make_view<double, 1>(matrix.value());

        ATLAS_ASSERT(outer_view(0) == 0);
        ATLAS_ASSERT(outer_view(1) == 1);
        ATLAS_ASSERT(outer_view(2) == 3);
        ATLAS_ASSERT(outer_view(3) == 3);
        ATLAS_ASSERT(outer_view(4) == 5);

        ATLAS_ASSERT(inner_view(0) == 1);
        ATLAS_ASSERT(inner_view(1) == 0);
        ATLAS_ASSERT(inner_view(2) == 1);
        ATLAS_ASSERT(inner_view(3) == 0);
        ATLAS_ASSERT(inner_view(4) == 2);

        ATLAS_ASSERT(values_view(0) == 2.);
        ATLAS_ASSERT(values_view(1) == 1.);
        ATLAS_ASSERT(values_view(2) == 3.);
        ATLAS_ASSERT(values_view(3) == 4.);
        ATLAS_ASSERT(values_view(4) == 5.);
    };

    const auto shuffle = [](const auto& begin, const auto& end) { std::shuffle(begin, end, std::mt19937{0}); };

    SECTION("Test sparse matrix construction from std::vector<TripletType>.") {
        auto [triplets, n_rows, n_cols] = make_test_triplets();
        // Triplets do not need to be sorted.
        test_triplets(make_sparse_matrix_storage_from_triplets(n_rows, n_cols, triplets, true));

        shuffle(triplets.begin(), triplets.end());
        // Triplets implicitly sorted.
        test_triplets(make_sparse_matrix_storage_from_triplets(n_rows, n_cols, triplets));


        // Should fail as triplets not sorted.
        shuffle(triplets.begin(), triplets.end());
        EXPECT_THROWS(make_sparse_matrix_storage_from_triplets(n_rows, n_cols, triplets, true));
    }

    SECTION("Test sparse matrix construction from TripletType[].") {
        auto [triplets, n_rows, n_cols] = make_test_triplets();
        test_triplets(make_sparse_matrix_storage_from_triplets(n_rows, n_cols, triplets.size(), triplets.data(), true));

        shuffle(triplets.begin(), triplets.end());
        test_triplets(make_sparse_matrix_storage_from_triplets(n_rows, n_cols, triplets.size(), triplets.data()));

        shuffle(triplets.begin(), triplets.end());
        EXPECT_THROWS(make_sparse_matrix_storage_from_triplets(n_rows, n_cols, triplets.size(), triplets.data(), true));
    }

    SECTION("Test sparse matrix construction from non-random-access-containers.") {
        auto [triplets, n_rows, n_cols] = make_test_triplets();

        const auto list_triplets = std::list<TripletType>{triplets.begin(), triplets.end()};

        shuffle(triplets.begin(), triplets.end());
        const auto shuffled_list_triplets = std::list<TripletType>{triplets.begin(), triplets.end()};

        // Note: std::set is sorted by construction.
        const auto set_triplets = std::set<TripletType>{triplets.begin(), triplets.end()};

        test_triplets(
            make_sparse_matrix_storage_from_triplets(n_rows, n_cols, list_triplets.begin(), list_triplets.end()));

        test_triplets(
            make_sparse_matrix_storage_from_triplets(n_rows, n_cols, set_triplets.begin(), set_triplets.end()));

        EXPECT_THROWS(make_sparse_matrix_storage_from_triplets(n_rows, n_cols, shuffled_list_triplets.begin(),
                                                               shuffled_list_triplets.end()));
    }

    SECTION("Test incorrect rows/cols throw") {
        auto [triplets, n_rows, n_cols] = make_test_triplets();
        shuffle(triplets.begin(), triplets.end());
        EXPECT_THROWS(make_sparse_matrix_storage_from_triplets(n_rows - 1, n_cols, triplets));
        EXPECT_THROWS(make_sparse_matrix_storage_from_triplets(n_rows, n_cols - 1, triplets));
    }
}

CASE("Test 'for each triplet' methods") {
    auto [triplets, n_rows, n_cols] = make_test_triplets();
    const auto ref_triplets         = triplets;
    const auto matrix_storage       = make_sparse_matrix_storage_from_triplets(n_rows, n_cols, triplets);
    const auto matrix_view          = make_host_view<double, int>(matrix_storage);

    SECTION("Test 'for each triplet' on matrix row.") {
        {
            auto ref_iter = ref_triplets.begin();
            for (std::size_t row = 0; row < matrix_view.rows(); ++row) {
                sparse_matrix_for_each_triplet(row, matrix_view, [&](int row, int col, double value) {
                    ATLAS_ASSERT(row == ref_iter->row());
                    ATLAS_ASSERT(col == ref_iter->col());
                    ATLAS_ASSERT(value == ref_iter->value());
                    ++ref_iter;
                });
            }
        }
    }
    SECTION("Test 'for each triplet' method on full matrix.") {
        {
            auto ref_iter = ref_triplets.begin();
            sparse_matrix_for_each_triplet(matrix_view, [&](int row, int col, double value) {
                ATLAS_ASSERT(row == ref_iter->row());
                ATLAS_ASSERT(col == ref_iter->col());
                ATLAS_ASSERT(value == ref_iter->value());
                ++ref_iter;
            });
        }
    }
}

}  // namespace atlas::test

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
