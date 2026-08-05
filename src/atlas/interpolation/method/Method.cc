/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <memory>
#include <vector>

#include "atlas/interpolation/method/Method.h"


#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/linalg/sparse.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

using namespace atlas::linalg;
namespace atlas {
namespace interpolation {

namespace {

template <typename Value>
void set_missing_values_rank1(Field& tgt, const std::vector<idx_t>& missing, const Value& missing_value) {
    auto tgt_v = array::make_view<Value, 1>(tgt);
    for (auto i : missing) {
        tgt_v(i) = missing_value;
    }
}

template <typename Value>
void set_missing_values_rank2(Field& tgt, const std::vector<idx_t>& missing, const Value& missing_value) {
    auto tgt_v     = array::make_view<Value, 2>(tgt);
    const idx_t Nj = tgt_v.shape(1);
    for (auto i : missing) {
        for (idx_t j = 0; j < Nj; ++j) {
            tgt_v(i, j) = missing_value;
        }
    }
}

template <typename Value>
void set_missing_values_rank3(Field& tgt, const std::vector<idx_t>& missing, const Value& missing_value) {
    auto tgt_v     = array::make_view<Value, 3>(tgt);
    const idx_t Nj = tgt_v.shape(1);
    const idx_t Nk = tgt_v.shape(2);
    for (auto i : missing) {
        for (idx_t j = 0; j < Nj; ++j) {
            for (idx_t k = 0; k < Nk; ++k) {
                tgt_v(i, j, k) = missing_value;
            }
        }
    }
}

template <typename Value>
void set_missing_values_T(Field& tgt, const std::vector<idx_t>& missing) {
    Value missing_value = tgt.metadata().get<Value>("missing_value");
    if (tgt.rank() == 1) {
        set_missing_values_rank1(tgt, missing, missing_value);
    }
    else if (tgt.rank() == 2) {
        set_missing_values_rank2(tgt, missing, missing_value);
    }
    else if (tgt.rank() == 3) {
        set_missing_values_rank3(tgt, missing, missing_value);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void set_missing_values(Field& tgt, const std::vector<idx_t>& missing) {
    if (missing.empty()) {
        return;
    }
    if (tgt.datatype().kind() == array::DataType::KIND_REAL64) {
        set_missing_values_T<double>(tgt, missing);
    }
    else if (tgt.datatype().kind() == array::DataType::KIND_REAL32) {
        set_missing_values_T<float>(tgt, missing);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

bool executesOnDevice(const sparse::Backend& backend) {
    if (backend.type() == "eckit_linalg") {
        return false;
    } else if (backend.type() == "openmp") {
        return false;
    } else if (backend.type() == "hicsparse") {
        return true;
    } else {
        ATLAS_NOTIMPLEMENTED;
    }
}

template<typename Value, int Rank>
atlas::array::ArrayView<Value, Rank> make_device_view_rw(atlas::Field& f) {
    f.syncDevice();
    return atlas::array::make_device_view<Value, Rank>(f);
}

template<typename Value, int Rank>
atlas::array::ArrayView<Value, Rank> make_device_view_w(atlas::Field& f) {
    if (not f.deviceAllocated()) {
        f.allocateDevice();
    }
    f.setDeviceNeedsUpdate(false); // as it will be written to here
    return atlas::array::make_device_view<Value, Rank>(f);
}

template<typename Value, int Rank>
atlas::array::ArrayView<const Value, Rank> make_device_view_r(const atlas::Field& f) {
    f.syncDevice();
    return atlas::array::make_device_view<const Value, Rank>(f);
}


template<typename Value, int Rank>
atlas::array::ArrayView<Value, Rank> make_host_view_rw(atlas::Field& f) {
    f.syncHost();
    return atlas::array::make_host_view<Value, Rank>(f);
}

template<typename Value, int Rank>
atlas::array::ArrayView<const Value, Rank> make_host_view_r(const atlas::Field& f) {
    f.syncHost();
    return atlas::array::make_host_view<const Value, Rank>(f);
}

template<typename Value, int Rank>
atlas::array::ArrayView<Value, Rank> make_host_view_w(atlas::Field& f) {
    f.setHostNeedsUpdate(false); // as it will be written to
    return atlas::array::make_host_view<Value, Rank>(f);
}

template<typename Value, typename Index>
atlas::linalg::SparseMatrixView<Value,Index> make_device_view_r(const atlas::linalg::SparseMatrixStorage& m) {
    if (m.deviceNeedsUpdate()) {
        ATLAS_TRACE("Copy interpolation matrix to device");
        m.updateDevice();
    }
    return make_device_view<Value, Index>(m);
}

template<typename Value, typename Index>
atlas::linalg::SparseMatrixView<Value,Index> make_host_view_r(const atlas::linalg::SparseMatrixStorage& m) {
    m.syncHost();
    return make_host_view<Value, Index>(m);
}

}  // anonymous namespace


template <typename Value>
void Method::interpolate_field_rank1(const Field& src, Field& tgt, const Matrix& W) const {
    auto backend = sparse::Backend{linalg_backend_};
    
    if (backend.type() == "hicsparse" && !std::is_same<eckit::linalg::Scalar, Value>::value) {
        ATLAS_NOTIMPLEMENTED; // hicsparse does not support mixed double-float
    }

    if (backend.type() == "eckit_linalg" && std::is_same<Value, float>::value) {
        // Switch to OpenMP as eckit_linalg does not support float
        backend = sparse::backend::openmp();
    }
    
    const auto on_device = executesOnDevice(backend);

    auto src_v = on_device ? make_device_view_r<Value, 1>(src) : make_host_view_r<Value, 1>(src);
    auto tgt_v = on_device ? make_device_view_w<Value, 1>(tgt) : make_host_view_w<Value, 1>(tgt);

    if (nonLinear_(src)) {
        eckit::linalg::SparseMatrix W_copy = atlas::linalg::make_eckit_sparse_matrix(W);
        nonLinear_->execute(W_copy, src);
        auto W_nl = make_sparse_matrix_storage(std::move(W_copy));
        auto W_nl_v = on_device ? make_device_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W_nl)
                                : make_host_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W_nl);
        sparse_matrix_multiply(W_nl_v, src_v, tgt_v, backend);
    }
    else {
        auto W_v = on_device ? make_device_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W)
                             : make_host_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W);
        sparse_matrix_multiply(W_v, src_v, tgt_v, backend);
    }

    on_device ? tgt.setHostNeedsUpdate(true) : tgt.setDeviceNeedsUpdate(true);
}


template <typename Value>
void Method::interpolate_field_rank2(const Field& src, Field& tgt, const Matrix& W) const {
    auto backend = sparse::Backend{linalg_backend_};
    
    if (backend.type() == "hicsparse" && !std::is_same<eckit::linalg::Scalar, Value>::value) {
        ATLAS_NOTIMPLEMENTED; // hicsparse does not support mixed double-float
    }

    if (backend.type() == "eckit_linalg") {
        // Switch to OpenMP as eckit_linalg does not support this layout
        backend = sparse::backend::openmp();
    }

    if (nonLinear_(src)) {
        // We cannot apply the same matrix to full columns as e.g. missing values could be present in only certain parts.
        //
        // Both the non-linear missing-value redistribution (nonlinear/Missing.cc) and the subsequent sparse matmul
        // (interpolate_field_rank1) only ever read the source at the columns referenced by the weights `W` (the stencil
        // support). We therefore only need to populate those referenced source rows in the per-level temporary, not the
        // entire source column. This turns the per-call copy from O(source_size * levels) into O(nnz * levels) ==
        // O(target_size * stencil * levels); for target sizes much smaller than the source (e.g. observation operators)
        // this eliminates the dominant, obs-independent apply cost.
        //
        // NOTE (possible future optimisation - caching W_nl): when the source missing-value *pattern* is static across
        // successive applies (e.g. a fixed land-sea / bathymetry mask), the per-level corrected matrices W_nl[lev]
        // produced by nonLinear_->execute() are constant and could be precomputed once at setup and cached. Apply could
        // then take the plain rank-2 sparse_matrix_multiply path (the `else` branch below) with no per-call
        // redistribution, no per-level matrix copy and no source copy at all - a further ~2-3x on the residual
        // obs-bound cost. That would need an opt-in flag (e.g. `non_linear_cache`) since atlas cannot assume every field
        // pushed through an interpolator shares one mask, plus storage for the per-level matrices (mitigable by
        // de-duplicating levels that share a mask). Left out here to keep this change assumption-free and local.
        //
        // NOTE (possible future optimisation - matmul source-read locality via compaction at construction): the gather
        // below only removes the *copy* cost; the subsequent matmul still reads `src_slice` scattered across a full
        // source-sized array (src_slice is allocated at src.shape(0) and W's column indices index the full source), so
        // it touches few useful entries per cache line. When target_size << source_size (observation operators), the
        // referenced-column set is small, so at setup one could build a compacted matrix W' whose columns are remapped
        // to a dense [0, |referenced|) space and cache it alongside the referenced-row list. Apply would then gather the
        // source into a small *dense* buffer of size |referenced| and matmul over W', shrinking the matmul working set
        // from source_size to |referenced| (fits in cache, full cache-line use). This needs nothing new at construction
        // (it is a pure structural remap of W, which is fully assembled in do_setup); it is static and cacheable.

        // Allocate temporary rank-1 fields corresponding to one horizontal level
        auto src_slice = Field("s", array::make_datatype<Value>(), {src.shape(0)});
        auto tgt_slice = Field("t", array::make_datatype<Value>(), {tgt.shape(0)});

        // Copy metadata to the source rank-1 field
        src_slice.metadata() = src.metadata();

        auto src_v = make_host_view_r<Value,2>(src);
        auto tgt_v = make_host_view_w<Value,2>(tgt);

        auto src_slice_v = array::make_host_view<Value, 1>(src_slice);
        auto tgt_slice_v = array::make_host_view<Value, 1>(tgt_slice);

        // For observation-operator-like interpolations (target << source) most source rows are unreferenced, so
        // gathering only the referenced rows per level avoids copying the full source column. For regridding-like
        // interpolations (target ~ source, nearly every source row referenced) the referenced set approaches the full
        // source and the extra sort/dedup below would be pure overhead. Use nnz vs source size as a cheap proxy (an
        // average of >= 1 nonzero per source row means most of the source is referenced) to fall back to the plain
        // full-column copy in that regime, keeping the regridding path free of the O(nnz log nnz) set construction.
        auto W_v = make_host_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W);
        const bool gather_referenced = static_cast<idx_t>(W_v.nnz()) < src.shape(0);

        // Unique set of source rows (= matrix columns) referenced by the weights; only built on the gather path.
        std::vector<eckit::linalg::Index> referenced_rows;
        if (gather_referenced) {
            referenced_rows.assign(W_v.inner(), W_v.inner() + W_v.nnz());
            std::sort(referenced_rows.begin(), referenced_rows.end());
            referenced_rows.erase(std::unique(referenced_rows.begin(), referenced_rows.end()), referenced_rows.end());
        }

        for (idx_t lev = 0; lev < src_v.shape(1); ++lev) {
            // Copy this level to the temporary rank-1 field: only the referenced source rows when gathering,
            // otherwise the full source column.
            if (gather_referenced) {
                for (const auto row : referenced_rows) {
                    src_slice_v(row) = src_v(row, lev);
                }
            }
            else {
                for (idx_t i = 0; i < src.shape(0); ++i) {
                    src_slice_v(i) = src_v(i, lev);
                }
            }

            // Interpolate between rank-1 fields
            interpolate_field_rank1<Value>(src_slice, tgt_slice, W);

            // Copy rank-1 field to this level in the rank-2 field
            tgt_slice.syncHost();
            for (idx_t i = 0; i < tgt.shape(0); ++i) {
                tgt_v(i, lev) = tgt_slice_v(i);
            }
        }
        tgt.setDeviceNeedsUpdate(true);
        tgt.setHostNeedsUpdate(false);
    }
    else {
        const auto on_device = executesOnDevice(backend);

        auto W_v = on_device ? make_device_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W)
                             : make_host_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W);
        auto src_dv = on_device ? make_device_view_r<Value, 2>(src) : make_host_view_r<Value, 2>(src);
        auto tgt_dv = on_device ? make_device_view_w<Value, 2>(tgt) : make_host_view_w<Value, 2>(tgt);
        
        sparse_matrix_multiply(W_v, src_dv, tgt_dv, backend);

        on_device ? tgt.setHostNeedsUpdate(true) : tgt.setDeviceNeedsUpdate(true);
    }
}


template <typename Value>
void Method::interpolate_field_rank3(const Field& src, Field& tgt, const Matrix& W) const {
    sparse::Backend backend{linalg_backend_};

    if (backend.type() == "hicsparse") {
        ATLAS_NOTIMPLEMENTED; // hicsparse does not support rank-3 fields
    }

    if (backend.type() == "eckit_linalg") {
        // Switch to OpenMP as eckit_linalg does not support rank-3 fields
        backend = sparse::backend::openmp();
    }

    auto W_v = make_host_view<eckit::linalg::Scalar, eckit::linalg::Index>(W);
    auto src_v = make_host_view_r<Value, 3>(src);
    auto tgt_v = make_host_view_w<Value, 3>(tgt);
    if (not W.empty() && nonLinear_(src)) {
        ATLAS_ASSERT(false, "nonLinear interpolation not supported for rank-3 fields.");
    }
    sparse_matrix_multiply(W_v, src_v, tgt_v, backend);

    tgt.setHostNeedsUpdate(false);
    tgt.setDeviceNeedsUpdate(true);
}

template <typename Value>
void Method::adjoint_interpolate_field_rank1(Field& src, const Field& tgt, const Matrix& W) const {
    auto backend = sparse::Backend{linalg_backend_};
    
    if (backend.type() == "hicsparse" && !std::is_same<eckit::linalg::Scalar, Value>::value) {
        ATLAS_NOTIMPLEMENTED; // hicsparse does not support mixed double-float
    }

    if (backend.type() == "eckit_linalg" && std::is_same<Value, float>::value) {
        // Switch to OpenMP as eckit_linalg does not support float
        backend = sparse::backend::openmp();
    }

    const auto on_device = executesOnDevice(backend);

    auto src_v = on_device ? make_device_view_rw<Value, 1>(src) : make_host_view_rw<Value, 1>(src);
    auto tgt_v = on_device ? make_device_view_r<Value, 1>(tgt) : make_host_view_r<Value, 1>(tgt);
    auto W_v = on_device ? make_device_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W)
                         : make_host_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W);

    sparse_matrix_multiply_add(W_v, tgt_v, src_v, backend);

    on_device ? src.setHostNeedsUpdate(true) : src.setDeviceNeedsUpdate(true);
}

template <typename Value>
void Method::adjoint_interpolate_field_rank2(Field& src, const Field& tgt, const Matrix& W) const {
    auto backend = sparse::Backend{linalg_backend_};
    
    if (backend.type() == "hicsparse" && !std::is_same<eckit::linalg::Scalar, Value>::value) {
        ATLAS_NOTIMPLEMENTED; // hicsparse does not support mixed double-float
    }

    if (backend.type() == "eckit_linalg") {
        // Switch to OpenMP as eckit_linalg does not support this layout
        backend = sparse::backend::openmp();
    }

    const auto on_device = executesOnDevice(backend);

    auto src_v = on_device ? make_device_view_rw<Value, 2>(src) : make_host_view_rw<Value, 2>(src);
    auto tgt_v = on_device ? make_device_view_r<Value, 2>(tgt) : make_host_view_r<Value, 2>(tgt);
    auto W_v = on_device ? make_device_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W)
                         : make_host_view_r<eckit::linalg::Scalar, eckit::linalg::Index>(W);

    sparse_matrix_multiply_add(W_v, tgt_v, src_v, backend);

    on_device ? src.setHostNeedsUpdate(true) : src.setDeviceNeedsUpdate(true);
}

template <typename Value>
void Method::adjoint_interpolate_field_rank3(Field& src, const Field& tgt, const Matrix& W) const {
    sparse::Backend backend{linalg_backend_};

    if (backend.type() == "hicsparse") {
        ATLAS_NOTIMPLEMENTED; // hicsparse does not support rank-3 fields
    }

    if (backend.type() == "eckit_linalg") {
        // Switch to OpenMP as eckit_linalg does not support rank-3 fields
        backend = sparse::backend::openmp();
    }

    auto src_v = make_host_view_rw<Value, 3>(src);
    auto tgt_v = make_host_view_r<Value, 3>(tgt);
    auto W_v = make_host_view_r<eckit::linalg::Scalar,eckit::linalg::Index>(W);

    sparse_matrix_multiply_add(W_v, tgt_v, src_v, backend);

    src.setDeviceNeedsUpdate(true);
}

void Method::check_compatibility(const Field& src, const Field& tgt, const Matrix& W) const {
    ATLAS_ASSERT(src.datatype() == tgt.datatype());
    ATLAS_ASSERT(src.rank() == tgt.rank());
    ATLAS_ASSERT(src.levels() == tgt.levels());
    ATLAS_ASSERT(src.variables() == tgt.variables());

    ATLAS_ASSERT(!W.empty());
    ATLAS_ASSERT(tgt.shape(0) >= static_cast<idx_t>(W.rows()));
    ATLAS_ASSERT(src.shape(0) >= static_cast<idx_t>(W.cols()));
}

template <typename Value>
void Method::interpolate_field(const Field& src, Field& tgt, const Matrix& W) const {
    // do nothing if there are no observations to interpolate (W will be NULL
    // and would fail the compatibility check)
    if (tgt.shape(0) == 0) {
        return;
    }
    check_compatibility(src, tgt, W);

    if (src.rank() == 1) {
        interpolate_field_rank1<Value>(src, tgt, W);
    }
    else if (src.rank() == 2) {
        interpolate_field_rank2<Value>(src, tgt, W);
    }
    else if (src.rank() == 3) {
        interpolate_field_rank3<Value>(src, tgt, W);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

template <typename Value>
void Method::adjoint_interpolate_field(Field& src, const Field& tgt, const Matrix& W) const {
    // do nothing if there are no observations to interpolate (W will be NULL
    // and would fail the compatibility check)
    if (tgt.shape(0) == 0) {
        return;
    }
    check_compatibility(tgt, src, W);

    if (src.rank() == 1) {
        adjoint_interpolate_field_rank1<Value>(src, tgt, W);
    }
    else if (src.rank() == 2) {
        adjoint_interpolate_field_rank2<Value>(src, tgt, W);
    }
    else if (src.rank() == 3) {
        adjoint_interpolate_field_rank3<Value>(src, tgt, W);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

Method::Method(const Method::Config& config) {
    config.get("sparse_matrix_multiply", linalg_backend_);  // empty is allowed -> sparse::current_backend()

    std::string non_linear;
    if (config.get("non_linear", non_linear)) {
        nonLinear_ = NonLinear(non_linear, config);
    }

    config.get("adjoint", adjoint_ = false);
}

void Method::setup(const FunctionSpace& source, const FunctionSpace& target) {
    ATLAS_TRACE("atlas::interpolation::method::Method::setup(FunctionSpace, FunctionSpace)");
    this->do_setup(source, target);
    if (adjoint_) {
        adjoint_matrix();
    }
}

const Method::Matrix& Method::adjoint_matrix() const {
    if (not matrix_transpose_) {
        if (target().size() == 0) {
            matrix_transpose_ = std::make_unique<Matrix>(); // Empty matrix
        }
        else {
            ATLAS_ASSERT(matrix_);
            eckit::linalg::SparseMatrix matrix_copy = make_eckit_sparse_matrix(*matrix_); // Makes a copy!
            matrix_copy.transpose(); // transpose the copy in place
            matrix_transpose_ = std::make_unique<Matrix>(linalg::make_sparse_matrix_storage(std::move(matrix_copy))); // Move the copy into storage
        }
    }
    return *matrix_transpose_;
}

void Method::setup(const Grid& source, const Grid& target) {
    ATLAS_TRACE("atlas::interpolation::method::Method::setup(Grid, Grid)");
    this->do_setup(source, target, Cache());
}

void Method::setup(const FunctionSpace& source, const Field& target) {
    ATLAS_TRACE("atlas::interpolation::method::Method::setup(FunctionSpace, Field)");
    this->do_setup(source, target);
}

void Method::setup(const FunctionSpace& source, const FieldSet& target) {
    ATLAS_TRACE("atlas::interpolation::method::Method::setup(FunctionSpace, FieldSet)");
    this->do_setup(source, target);
}

void Method::setup(const Grid& source, const Grid& target, const Cache& cache) {
    ATLAS_TRACE("atlas::interpolation::method::Method::setup(Grid, Grid, Cache)");
    this->do_setup(source, target, cache);
}

void Method::setup(const FunctionSpace& source, const FunctionSpace& target, const Cache& cache) {
    ATLAS_TRACE("atlas::interpolation::method::Method::setup(FunctionSpace, FunctionSpace, Cache)");
    this->do_setup(source, target, cache);
}

Method::Metadata Method::execute(const FieldSet& source, FieldSet& target) const {
    ATLAS_TRACE("atlas::interpolation::method::Method::execute(FieldSet, FieldSet)");
    Metadata metadata;
    this->do_execute(source, target, metadata);
    return metadata;
}

Method::Metadata Method::execute(const Field& source, Field& target) const {
    ATLAS_TRACE("atlas::interpolation::method::Method::execute(Field, Field)");
    Metadata metadata;
    this->do_execute(source, target, metadata);
    return metadata;
}

Method::Metadata Method::execute_adjoint(FieldSet& source, const FieldSet& target) const {
    ATLAS_TRACE("atlas::interpolation::method::Method::execute_adjoint(FieldSet, FieldSet)");
    Metadata metadata;
    this->do_execute_adjoint(source, target, metadata);
    return metadata;
}

Method::Metadata Method::execute_adjoint(Field& source, const Field& target) const {
    ATLAS_TRACE("atlas::interpolation::method::Method::execute_adjoint(Field, Field)");
    Metadata metadata;
    this->do_execute_adjoint(source, target, metadata);
    return metadata;
}

void Method::do_setup(const FunctionSpace& /*source*/, const Field& /*target*/) {
    ATLAS_NOTIMPLEMENTED;
}

void Method::do_setup(const FunctionSpace& /*source*/, const FieldSet& /*target*/) {
    ATLAS_NOTIMPLEMENTED;
}

void Method::do_execute(const FieldSet& fieldsSource, FieldSet& fieldsTarget, Metadata& metadata) const {
    ATLAS_TRACE("atlas::interpolation::method::Method::do_execute()");

    const idx_t N = fieldsSource.size();
    ATLAS_ASSERT(N == fieldsTarget.size());

    for (idx_t i = 0; i < fieldsSource.size(); ++i) {
        Method::do_execute(fieldsSource[i], fieldsTarget[i], metadata);
    }
}

void Method::do_execute(const Field& src, Field& tgt, Metadata&) const {
    ATLAS_TRACE("atlas::interpolation::method::Method::do_execute()");

     if (src.hostNeedsUpdate() && src.deviceNeedsUpdate()) {
            throw_AssertionFailed("Inconsistent memory state flags for field "+src.name()+"  (hostNeedsUpdate && deviceNeedsUpdate) - we will not be able to "
                                  "determine which memory space to perform the halo exchange on",
                                  Here());
    }

    sparse::Backend backend{linalg_backend_};

    const bool on_device = executesOnDevice(backend);

    haloExchange(src, on_device);
    
    if( matrix_ ) { // (matrix == nullptr) when a partition is empty
        if (src.datatype().kind() == array::DataType::KIND_REAL64) {
            interpolate_field<double>(src, tgt, *matrix_);
        }
        else if (src.datatype().kind() == array::DataType::KIND_REAL32) {
            interpolate_field<float>(src, tgt, *matrix_);
        }
        else {
            ATLAS_NOTIMPLEMENTED;
        }
    }

    // carry over missing value metadata
    if (not tgt.metadata().has("missing_value")) {
        field::MissingValue mv_src(src);
        if (mv_src) {
            mv_src.metadata(tgt);
            ATLAS_ASSERT(field::MissingValue(tgt));
        }
        else if (not missing_.empty()) {
            if (not tgt.metadata().has("missing_value")) {
                tgt.metadata().set("missing_value", 9999.);
            }
            tgt.metadata().set("missing_value_type", "equals");
        }
    }

    // set missing values on host
    if (not missing_.empty()) {
        if (on_device) {
            tgt.updateHost();
        }
        set_missing_values(tgt, missing_);
        tgt.set_dirty();
        tgt.setDeviceNeedsUpdate(true);
    }
    else {
        tgt.set_dirty();
        on_device ? tgt.setHostNeedsUpdate(true) : tgt.setDeviceNeedsUpdate(true);
    }
}

void Method::do_execute_adjoint(FieldSet& fieldsSource, const FieldSet& fieldsTarget, Metadata& metadata) const {
    ATLAS_TRACE("atlas::interpolation::method::Method::do_execute_adjoint()");

    const idx_t N = fieldsSource.size();
    ATLAS_ASSERT(N == fieldsTarget.size());

    for (idx_t i = 0; i < fieldsSource.size(); ++i) {
        Method::do_execute_adjoint(fieldsSource[i], fieldsTarget[i], metadata);
    }
}

void Method::do_execute_adjoint(Field& src, const Field& tgt, Metadata&) const {
    ATLAS_TRACE("atlas::interpolation::method::Method::do_execute_adjoint()");

    if (src.hostNeedsUpdate() && src.deviceNeedsUpdate()) {
            throw_AssertionFailed("Inconsistent memory state flags - we will not be able to "
                                  "determine which memory space to perform the adjoint halo exchange on",
                                  Here());
    }

    sparse::Backend backend{linalg_backend_};

    if (nonLinear_(src)) {
        throw_NotImplemented("Adjoint interpolation only works for interpolation schemes that are linear", Here());
    }

    if (not missing_.empty()) {
        throw_NotImplemented("Adjoint Interpolation does not work for fields that have missing data. ", Here());
    }

    if (src.datatype().kind() == array::DataType::KIND_REAL64) {
        adjoint_interpolate_field<double>(src, tgt, adjoint_matrix());
    }
    else if (src.datatype().kind() == array::DataType::KIND_REAL32) {
        adjoint_interpolate_field<float>(src, tgt, adjoint_matrix());
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }

    src.set_dirty();

    const bool on_device = executesOnDevice(backend);

    adjointHaloExchange(src, on_device);
    
    on_device ? src.setHostNeedsUpdate(true) : src.setDeviceNeedsUpdate(true);
}


void Method::normalise(Triplets& triplets) {
    // sum all calculated weights for normalisation
    double sum = 0.0;

    for (size_t j = 0; j < triplets.size(); ++j) {
        sum += triplets[j].value();
    }

    // now normalise all weights according to the total
    const double invSum = 1.0 / sum;
    for (size_t j = 0; j < triplets.size(); ++j) {
        triplets[j].value() *= invSum;
    }
}

void Method::haloExchange(const FieldSet& fields, bool on_device) const {
    for (auto& field : fields) {
        haloExchange(field, on_device);
    }
}
void Method::haloExchange(const Field& field, bool on_device) const {
    if (field.dirty() && allow_halo_exchange_) {
        ATLAS_TRACE("haloExchange");
        if (on_device) {
            if (field.deviceNeedsUpdate()) {
                // Prefer halo exchange on host and copy to device
                source().haloExchange(field); // on host
                field.syncDevice();
            }
            else {
                source().haloExchange(field, on_device); // on device
            }
        }
        else {
            source().haloExchange(field); // on host
        }
    }
}

void Method::adjointHaloExchange(const FieldSet& fields, bool on_device) const {
    for (auto& field : fields) {
        adjointHaloExchange(field, on_device);
    }
}
void Method::adjointHaloExchange(const Field& field, bool on_device) const {
    if (field.dirty() && allow_halo_exchange_) {
        ATLAS_TRACE("adjointHaloExchange");
        if (on_device) {
            if (field.deviceNeedsUpdate()) {
                // Prefer halo exchange on host and copy to device
                source().adjointHaloExchange(field); // on host
                field.syncDevice();
            }
            else {
                source().adjointHaloExchange(field, on_device); // on device
            }
        }
        else {
            source().adjointHaloExchange(field); // on host
        }
    }
}

interpolation::Cache Method::createCache() const {
    return matrix_cache_;
}


}  // namespace interpolation
}  // namespace atlas
