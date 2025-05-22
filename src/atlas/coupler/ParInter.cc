/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/coupler/ParInter.h"


namespace atlas {


ParInter::ParInter(const Matrix& gmat, const FunctionSpace src_fs, const FunctionSpace tgt_fs,
    const grid::Distribution& src_distribution, const grid::Distribution& tgt_distribution) : 
    gmat_(gmat), src_fs_(src_fs), tgt_fs_(tgt_fs), src_distribution_(src_distribution), tgt_distribution_(tgt_distribution) {
        ATLAS_TRACE("Parinter::setup");
        extract(tgt_fs, tgt_distribution, gmat_, rows_, cols_, gcols_, vals_);
        setup_collect(src_fs_, gcols_);

        std::size_t nr = tgt_fs.size();
        std::size_t nc = collect_size_;
        eckit::linalg::Index index_base = 0;
        bool is_sorted = false;
        lmat_ = linalg::make_sparse_matrix_storage_from_rows_columns_values(nr, nc, rows_, cols_, vals_, index_base, is_sorted);
        find_missing_rows();

        nonlinear_ = interpolation::NonLinear("missing-if-any-missing", util::Config());

}


void ParInter::print() const {
    for (int task = 0; task < mpi::comm().size(); ++task) {
        if (mpi::comm().rank() == task) {
            std::cout << "TASK " << task << std::endl;
            for (int i = 0; i < rows_.size(); ++i) {
                std::cout << "\t" << rows_[i] << ", " << cols_[i] << ", " << vals_[i];
                if (i < gcols_.size()) {
                    std::cout << ", " << gcols_[i];
                }
            }
        }
        mpi::comm().barrier();
    }
}


void ParInter::execute(const Field& src, Field& tgt) {
    ATLAS_TRACE("Parinter::execute");
    auto collect_shape = src.shape();
    collect_shape[0] = collect_size_;
    std::unique_ptr<array::Array> collect_src(array::Array::create(src.datatype(), collect_shape));

    collect_.execute<double, 1>(const_cast<Field&>(src), *collect_src);

    {
        auto collect_src_view = array::make_view<double,1>(*collect_src);
        auto tgt_view = array::make_view<double,1>(tgt);
        auto matrix = linalg::make_host_view<eckit::linalg::Scalar, eckit::linalg::Index>(lmat_);

        ATLAS_ASSERT(tgt_view.shape(0) == matrix.rows());

        ATLAS_TRACE_SCOPE("sparse_matrix_multiply") {
            // linalg::sparse::current_backend("eckit_linalg");
            if (nonlinear_(src)) {
                eckit::linalg::SparseMatrix matrix_nl = atlas::linalg::make_eckit_sparse_matrix(matrix);
                ATLAS_ASSERT(matrix_nl.cols() == collect_src->size());
                ATLAS_ASSERT(matrix_nl.cols() == collect_src->size());
                nonlinear_->execute(matrix_nl, src, *collect_src);
                linalg::sparse_matrix_multiply(matrix_nl, collect_src_view, tgt_view);
            }
            else {
                linalg::sparse_matrix_multiply(matrix, collect_src_view, tgt_view);
            }
        }

        if (not tgt.metadata().has("missing_value")) {
            field::MissingValue mv_src(src);
            if (mv_src) {
                ATLAS_DEBUG();
                mv_src.metadata(tgt);
                ATLAS_ASSERT(field::MissingValue(tgt));
            }
            else if (not missing_rows_.empty()) {
                ATLAS_DEBUG();
                if (not tgt.metadata().has("missing_value")) {
                    tgt.metadata().set("missing_value", 9999.);
                }
                tgt.metadata().set("missing_value_type", "approximately-equals");
                tgt.metadata().set("missing_value_epsilon", 1.);
            }
            tgt.metadata().set("missing_value", 9999.);
            tgt.metadata().set("missing_value_type", "approximately-equals");
            tgt.metadata().set("missing_value_epsilon", 1.);

        }
        
        ATLAS_TRACE_SCOPE("mask") {
            ATLAS_DEBUG();
            if (tgt.metadata().has("missing_value")) {
                double missing_value = tgt.metadata().get<double>("missing_value");
                for (idx_t r : missing_rows_) {
                    tgt_view(r) = missing_value; 
                }
            }
        }
    }
}


void ParInter::find_missing_rows() {
    auto matrix = linalg::make_host_view<eckit::linalg::Scalar, eckit::linalg::Index>(lmat_);
    for(std::size_t r = 0; r < matrix.rows(); ++r) {
        int cols = matrix.outer()[r + 1] - matrix.outer()[r];
        if (cols == 0) {
            missing_rows_.emplace_back(r);
        }
    }
}


} // namespace atlas