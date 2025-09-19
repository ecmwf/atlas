/*
 * (C) Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <random>
#include "atlas/linalg/sparse/SparseMatrixStorage.h"
#include "atlas/linalg/sparse/SparseMatrixTriplet.h"
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
};

CASE("Test sparse matrix construction from triplets.") {
    SECTION("Test sparse matrix construction from triplets.") {
        auto [triplets, n_rows, n_cols] = make_test_triplets();
        // Triplets do not need to be sorted.
        std::shuffle(triplets.begin(), triplets.end(), std::default_random_engine(0));

        const auto matrix = make_sparse_matrix_storage_from_triplets(n_rows, n_cols, std::move(triplets));

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
    }

    SECTION("Test 'for each triplet' method on sparse matrix.") {
        auto [triplets, n_rows, n_cols] = make_test_triplets();
        const auto ref_triplets         = triplets;
        const auto matrix_storage       = make_sparse_matrix_storage_from_triplets(n_rows, n_cols, std::move(triplets));
        const auto matrix_view          = make_host_view<double, int>(matrix_storage);

        {
            // Row-wise iteration.
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
        {
            // Full matrix iteration.
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
