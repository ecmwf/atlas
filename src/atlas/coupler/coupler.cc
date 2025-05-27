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
#include "atlas/coupler/ModelCoupler.h"
#include "atlas/coupler/ParInter.h"
#include "atlas/parallel/mpi/mpi.h"


void setup_oneway_remap(int model_1, int model_2) {
    using namespace atlas;
    using namespace atlas::coupler;
    using Matrix      = atlas::linalg::SparseMatrixStorage;

    // create the model1 & model2 joint communicator
    mpi::comm("world").barrier();
    auto m1_ranks = coupler_.global_ranks(10);
    auto m2_ranks = coupler_.global_ranks(21);
    auto in_m1_m2_ranks = [&m1_ranks, &m2_ranks](int rank) {
        bool in_m1 = (std::find(m1_ranks.begin(), m1_ranks.end(), rank) != m1_ranks.end());
        if (in_m1) {
            return true;
        }
        bool in_m2 = (std::find(m2_ranks.begin(), m2_ranks.end(), rank) != m2_ranks.end());
        return in_m2;
    };
    auto& comm12 = eckit::mpi::comm("world").split(in_m1_m2_ranks(mpi::comm("world").rank()), "comm12");
    mpi::comm("world").barrier();

    bool is_model_1 = (mpi::comm().name() == coupler_.comm_name(model_1));
    bool is_model_2 = (mpi::comm().name() == coupler_.comm_name(model_2));
    if (! is_model_1 && ! is_model_2) {
        return;
    }

    // exchange models
    std::vector<int> models;
    models.resize(comm12.size());
    models[comm12.rank()] = coupler_.this_model();
    comm12.allGather(coupler_.this_model(), models.begin(), models.end());

    // exchange other procs' comm12-ranks
    std::unordered_map<int, std::vector<int>> model_comm12_ranks;
    {
        std::vector<int> comm12_sizes(comm12.size());
        int my_size = mpi::comm().size();
        comm12.allGather(my_size, comm12_sizes.begin(), comm12_sizes.end());

        std::vector<int> my_comm12_ranks(my_size);
        int my_rank = comm12.rank();
        mpi::comm().allGather(my_rank, my_comm12_ranks.begin(), my_comm12_ranks.end());

        eckit::mpi::Buffer<int> buf_ranks(comm12.size());
        comm12.allGatherv(my_comm12_ranks.begin(), my_comm12_ranks.end(), buf_ranks);

        std::vector<int> comm12_models;
        for (int i = 0; i < comm12.size(); i++) {
            if (std::find(comm12_models.begin(), comm12_models.end(), models[i]) == comm12_models.end()) {
                comm12_models.emplace_back(models[i]);
                model_comm12_ranks[models[i]].reserve(buf_ranks.counts[i]);
                std::copy_n(&(buf_ranks.buffer[buf_ranks.displs[i]]), buf_ranks.counts[i],
                        std::back_inserter(model_comm12_ranks[models[i]]));
            }
        }
    }

    // everything, except the remote index, is set up on model_2 for remapping
    // global matrix is read in on model_2
    atlas::Grid grid1;
    atlas::Grid grid2; // later, no need to recreate grid2, as we are on the tasks which own grid2
    Matrix gmatrix;
    Matrix lmatrix;
    auto collect = coupler_.collect_map(model_1, model_2);

    if (is_model_2) {
        if (atlas::mpi::comm().rank() == 0) {
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
            bool is_sorted = true; // !!! comming from Atlas itself, the tripplets are sorted
            lmatrix = linalg::make_sparse_matrix_storage_from_rows_columns_values(nr, nc, l_rows, l_cols, l_vals, index_base, is_sorted);
        }

        grid1 = atlas::Grid(coupler_.grid_name(model_1));
        grid::Distribution dist_m1(grid1, grid1.partitioner()
            | util::Config("nb_partitions", coupler_.global_ranks(model_1).size())
            | util::Config("mpi_comm", "21"));

        std::size_t collect_size = l_gcols.size();
        std::vector<gidx_t> collect_gidx(collect_size);
        std::vector<idx_t>  collect_ridx(collect_size);
        std::vector<int> collect_partition(collect_size);
        idx_t ridx_base = 0;

        for (size_t i = 0; i < collect_gidx.size(); ++i) {
            collect_gidx[i] = l_gcols[i] + 1;
        }

        util::locate_partition(dist_m1, collect_gidx.size(), collect_gidx.data(), collect_partition.data());

        for (int i = 0; i < collect_partition.size(); ++i) {
            auto p = collect_partition[i]; // comm1-rank
            collect_partition[i] = model_comm12_ranks.at(model_1)[p]; // comm12-rank
        }

        util::locate_remote_index(comm12, 0, nullptr, nullptr,
                                    collect_size, collect_gidx.data(), collect_partition.data(), collect_ridx.data(), ridx_base);

        collect->setup(comm12.name(), collect_size, collect_partition.data(), collect_ridx.data(), ridx_base);
    }

    if (is_model_1) {
        grid1 = atlas::Grid(coupler_.grid_name(model_1));
        functionspace::StructuredColumns fspace_m1(grid1);
        auto glb_idx_v = array::make_view<gidx_t, 1>(fspace_m1.global_index());
        auto ghost_v   = array::make_view<int, 1>(fspace_m1.ghost());
        util::locate_remote_index(comm12, glb_idx_v.size(), glb_idx_v.data(), ghost_v.data(),
                                    0, nullptr, nullptr, nullptr, 0);
        collect->setup(comm12.name(), 0, nullptr, nullptr, 0);
    }
    comm12.barrier();
}


namespace atlas::coupler {


    void setup(int model_id, Grid grid) {
        eckit::mpi::comm().split(model_id, coupler_.comm_name(model_id));
        eckit::mpi::setCommDefault(coupler_.comm_name(model_id));
        coupler_.register_model(model_id, grid);
        coupler_.print_models();

        // HACK: we assume only m1 sends field 'f1' to m2
        setup_oneway_remap(10, 21);
    }


    void finalise() {
        coupler_.finalise();
    }


    void put(Field f, int model_2, int tstep) {
        auto& collect = *(coupler_.collect_map(coupler_.this_model(), model_2));
        collect.send<double, 1>(f);
    }


    void get(Field f, int model_1, int tstep) {
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
