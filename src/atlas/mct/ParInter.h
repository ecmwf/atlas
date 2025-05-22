/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/Iterator.h"
#include "atlas/interpolation.h"
#include "atlas/interpolation/AssembleGlobalMatrix.h"
#include "atlas/interpolation/NonLinear.h"
#include "atlas/io/atlas-io.h"
#include "atlas/linalg/sparse.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/linalg/sparse/SparseMatrixToTriplets.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/Collect.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/function/VortexRollup.h"
#include "atlas/util/function/MDPI_functions.h"
#include "atlas/util/Locate.h"
#include "atlas/field/MissingValue.h"


using atlas::functionspace::PointCloud;
using atlas::functionspace::NodeColumns;
using atlas::functionspace::StructuredColumns;
using atlas::util::Config;
using Matrix = atlas::linalg::SparseMatrixStorage;

namespace atlas {

class ParInter {
public:
    ParInter(const Matrix& gmat, const FunctionSpace src_fs, const FunctionSpace tgt_fs, 
        const grid::Distribution& src_distribution, const grid::Distribution& tgt_distribution);

    void print() const;

    void execute(const Field& src, Field& tgt);

    template <typename Vector>
    void setup_collect(FunctionSpace fs, const Vector& global_index);

    template <typename Value, typename Index> 
    static void extract(const FunctionSpace tgt_fspace, const grid::Distribution tgt_dist, const Matrix& gmat, std::vector<Index>& local_rows,
        std::vector<Index>& local_cols, std::vector<Index>& local_gcols, std::vector<Value>& local_vals, int mpi_root = 0);

    void find_missing_rows();

private:
    const Matrix& gmat_;
    const FunctionSpace src_fs_;
    const FunctionSpace tgt_fs_;
    const grid::Distribution src_distribution_;
    const grid::Distribution tgt_distribution_;
    parallel::Collect collect_;
    std::vector<eckit::linalg::Index> rows_;
    std::vector<eckit::linalg::Index> cols_;
    std::vector<eckit::linalg::Index> gcols_;
    std::vector<eckit::linalg::Scalar> vals_;
    Matrix lmat_;
    std::size_t collect_size_;

    interpolation::NonLinear nonlinear_;
    
    std::vector<idx_t> missing_rows_;

}; // class ParInter


template <typename Vector>
void ParInter::setup_collect(FunctionSpace fs, const Vector& global_index) {
        ATLAS_TRACE();
    collect_size_ = global_index.size();

    std::vector<gidx_t> collect_gidx(collect_size_);
    std::vector<int>    collect_partition(collect_size_);
    std::vector<idx_t>  collect_ridx(collect_size_);

    for (size_t i=0; i<global_index.size(); ++i) {
        collect_gidx[i] = global_index[i] + 1;
    }
    idx_t ridx_base = 0;

    ATLAS_ASSERT(fs);
    ATLAS_ASSERT(src_distribution_);
    ATLAS_DEBUG_VAR(collect_size_);
    ATLAS_DEBUG_VAR(src_distribution_.size());
    ATLAS_DEBUG_VAR(fs.size());
    util::locate(fs, src_distribution_, collect_gidx, collect_partition, collect_ridx, ridx_base);
    collect_.setup(collect_size_, collect_partition.data(), collect_ridx.data(), ridx_base);
}


template <typename Value, typename Index> 
void ParInter::extract(const FunctionSpace tgt_fspace, const grid::Distribution tgt_dist, const Matrix& gmat,
    std::vector<Index>& local_rows, std::vector<Index>& local_cols, std::vector<Index>& local_gcols,
    std::vector<Value>& local_vals, int mpi_root) {
    ATLAS_TRACE("ParInter::extract");

    // Field field_fs_part_glb = tgt_fs_.createField(tgt_fs_.partition(), option::global(mpi_root));
    // ATLAS_TRACE_SCOPE("gather partition") {
    //     tgt_fs_.gather(tgt_fs_.partition(), field_fs_part_glb);
    // }
    std::vector<Index> rows, cols;
    std::vector<Value> vals;
    interpolation::distribute_global_matrix_as_triplets(
        // array::make_view<int,1>(field_fs_part_glb).data(), 
        tgt_dist,
        atlas::linalg::make_host_view<Value, Index>(gmat), rows, cols, vals, mpi_root);

    std::unordered_map<gidx_t, idx_t> to_local_rows;
    ATLAS_TRACE_SCOPE("convert to local row indexing") {
        auto fs_gidx_exchanged = tgt_fspace.createField(tgt_fspace.global_index());
        fs_gidx_exchanged.array().copy(tgt_fspace.global_index());
        tgt_fspace.haloExchange(fs_gidx_exchanged);
        const auto fs_global_index = array::make_view<gidx_t, 1>(fs_gidx_exchanged);
        Field fs_ghost_field;
        if (functionspace::CellColumns(tgt_fspace)) {
            auto fs = functionspace::CellColumns(tgt_fspace);
            fs_ghost_field = fs.mesh().cells().halo();
        }
        else {
            fs_ghost_field = tgt_fspace.ghost();
        }
        const auto fs_ghost = array::make_view<int,1>(fs_ghost_field);

        for (idx_t r = 0; r < fs_global_index.size(); ++r) {
            auto gr = fs_global_index(r) - 1;
            if (fs_ghost(r) && to_local_rows.find(gr) != to_local_rows.end()) {
                continue;
            }
            to_local_rows[gr] = r;
        }

        std::unordered_map<gidx_t,idx_t> map_gcol_to_loc_col_idx;
        auto find_loc_col_idx = [&map_gcol_to_loc_col_idx](const gidx_t gcol) {
            auto it = map_gcol_to_loc_col_idx.find(gcol);
            if (it == map_gcol_to_loc_col_idx.end()) {
                return -1;
            }
            return it->second;
        };

        idx_t loc_col_idx = 0;
        auto new_loc_col_idx = [&](const gidx_t gcol) {
            idx_t lcol = loc_col_idx;
            map_gcol_to_loc_col_idx[gcol] = lcol;
            ++loc_col_idx;
            return lcol;
        };

        for(size_t idx = 0; idx < rows.size(); ++idx) {
            auto loc_row = to_local_rows[rows[idx]];
            auto glb_col = cols[idx];
            auto val = vals[idx];
            local_rows.emplace_back(loc_row);
            local_vals.emplace_back(val);
            auto found_loc_col_idx = find_loc_col_idx(glb_col);
            if (found_loc_col_idx >= 0) {
                local_cols.emplace_back(found_loc_col_idx);
            }
            else {
                local_cols.emplace_back(new_loc_col_idx(glb_col));
                local_gcols.emplace_back(glb_col);
            }
        }
    }
}


} // namespace atlas