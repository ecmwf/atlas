/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <iterator>
#include <vector>
#include <map>

#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/interpolation/AssembleGlobalMatrix.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/coupler/coupler.h"
#include "atlas/coupler/ParInter.h"
#include "atlas/parallel/mpi/mpi.h"

using Matrix      = atlas::linalg::SparseMatrixStorage;

namespace atlas::coupler {

    void set_default_comm_to_local(int model_id)
    {
        eckit::mpi::comm().split(model_id, coupler_.comm_name(model_id));
        eckit::mpi::setCommDefault(coupler_.comm_name(model_id));
    }


    void setup_oneway_remap(int model_1, int model_2) {
        bool is_model_1 = (mpi::comm().name() == coupler_.comm_name(model_1));
        bool is_model_2 = (mpi::comm().name() == coupler_.comm_name(model_2));
        if (! is_model_1 && ! is_model_2) {
            return;
        }

        // everything, except the remote index, is set up on model_2 for remapping
        // global matrix is read in on model_2

        atlas::Grid grid1;
        atlas::Grid grid2; // later, no need to recreate grid2, as we are on the tasks which own grid2

        Matrix gmatrix;
        Matrix lmatrix;
        auto collect = coupler_.collect_map(model_1, model_2);

        if (is_model_2) {
            if (mpi::comm().rank() == 0) {
                eckit::linalg::SparseMatrix eckit_matrix;
                eckit_matrix.load("mat.eckit");
                gmatrix = atlas::linalg::make_sparse_matrix_storage(std::move(eckit_matrix));
            }
            grid2 = atlas::Grid(coupler_.grid_name(model_2));
            grid::Distribution dist_m2(grid2, grid2.partitioner());
            functionspace::StructuredColumns fspace_m2(grid2);

            std::vector<Index> l_rows;
            std::vector<Index> l_cols;
            std::vector<Index> l_gcols;
            std::vector<Value> l_vals;
            ParInter::extract(fspace_m2, dist_m2, gmatrix, l_rows, l_cols, l_gcols, l_vals);

            auto& lmatrix = coupler_.remap(model_1, model_2);
            {
                std::size_t nr = fspace_m2.size();
                std::size_t nc = l_gcols.size();
                eckit::linalg::Index index_base = 0;
                bool is_sorted = true; // !!! check if tripplets are sorted here
                lmatrix = linalg::make_sparse_matrix_storage_from_rows_columns_values(nr, nc, l_rows, l_cols, l_vals, index_base, is_sorted);
            }

            grid1 = atlas::Grid(coupler_.grid_name(model_1));
            grid::Distribution dist_m1(grid1, grid1.partitioner() | util::Config("nb_partitions", coupler_.global_ranks(model_1).size()));

            std::size_t collect_size = l_gcols.size();
            std::vector<gidx_t> collect_gidx(collect_size);
            std::vector<idx_t>  collect_ridx(collect_size);
            std::vector<int> collect_partition;
            idx_t ridx_base = 0;

            for (size_t i=0; i < collect_gidx.size(); ++i) {
                collect_gidx[i] = l_gcols[i] + 1;
            }

            collect_partition.resize(collect_size);
            util::locate_partition(dist_m1, collect_gidx.size(), collect_gidx.data(), collect_partition.data());
            for (int i = 0; i < collect_partition.size(); ++i) {
                auto part = collect_partition[i]; // local rank
                collect_partition[i] = coupler_.global_ranks(model_1)[part]; // "world"-rank
            }
            util::locate_remote_index(mpi::comm("world"), 0, nullptr, nullptr,
                                      collect_size, collect_gidx.data(), collect_partition.data(), collect_ridx.data(), ridx_base);
            collect->setup("world", collect_size, collect_partition.data(), collect_ridx.data(), ridx_base);
            mpi::comm().barrier();
        }

        if (is_model_1) {
            grid1 = atlas::Grid(coupler_.grid_name(model_1));
            functionspace::StructuredColumns fspace_m1(grid1);
            auto glb_idx_v = array::make_view<gidx_t, 1>(fspace_m1.global_index());
            auto ghost_v   = array::make_view<int, 1>(fspace_m1.ghost());
            util::locate_remote_index(mpi::comm("world"), glb_idx_v.size(), glb_idx_v.data(), ghost_v.data(),
                                      0, nullptr, nullptr, nullptr, 0);
            collect->setup("world", 0, nullptr, nullptr, 0);
        }
    }


    void setup(int model_id, Grid grid) {
        set_default_comm_to_local(model_id);
        coupler_.register_model(model_id, grid);
        coupler_.print_models();

        // HACK: we assume only m1 sends field 'f1' to m2
        setup_oneway_remap(10, 21);
    }


    void finalise() {
        coupler_.finalise();
    }


    void put_field(Field f, int model_2, int tstep) {
        auto& collect = *(coupler_.collect_map(coupler_.this_model(), model_2));
        collect.send<double, 1>(f);
    }


    void get_field(Field f, int model_1, int tstep) {
        auto& collect = *(coupler_.collect_map(model_1, coupler_.this_model()));
        std::size_t size = collect.recv_size(); // !!! shape has to come from src_field
        std::unique_ptr<array::Array> collect_src{ array::Array::create(f.datatype(), array::make_shape(size)) };
        collect.recv<double, 1>(*collect_src);
        // interpolate to f
        const auto& lmatrix = coupler_.remap(model_1, coupler_.this_model());
        auto matrix = linalg::make_host_view<eckit::linalg::Scalar, eckit::linalg::Index>(lmatrix);
        const auto collect_src_view = array::make_view<double,1>(*collect_src);
        auto tgt_field_v = array::make_view<double, 1>(f);
        linalg::sparse_matrix_multiply(matrix, collect_src_view, tgt_field_v);
    }

} // namespace atlas::coupler
