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

#include <iterator>
#include <vector>
#include <map>

#include "atlas/grid.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/interpolation/AssembleGlobalMatrix.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/mct/ParInter.h"
#include "atlas/parallel/mpi/mpi.h"

using Matrix      = atlas::linalg::SparseMatrixStorage;
using Index = eckit::linalg::Index;
using Value = eckit::linalg::Scalar;

namespace atlas::mct {

    class ModelCoupler {
    public:
        using CollectMapType = std::map<std::string, std::unique_ptr<parallel::Collect>>;

        ModelCoupler() {}

        void register_model(int model_id, Grid grid);
        int this_model() const { return this_model_id_; }
        void print_models() const;

        const std::vector<int>& models() const { return models_; }
        const std::string comm_name(int model_id) const { return std::to_string(model_id); }

        const std::string& grid_name(int model_id) const { return model2grid_.at(model_id); }
        const std::vector<int>& global_ranks(int model_id) const { return model2ranks_.at(model_id); }

        void finalise();

        parallel::Collect* collect_map(int model_1, int model_2) {
            std::string key = std::to_string(model_1) + "_" + std::to_string(model_2);
            collect_map_.emplace(key, new parallel::Collect());
            return collect_map_.at(key);
        }

        // eckit::mpi::Request* recv_request(int model_1, int model_2, const Field f, int tstep) {
        //     std::string key = std::to_string(model_1) + "_" + std::to_string(model_2) 
        //         + "_f" + f.name() + "_t" + std::to_string(tstep);
        //     collect_map_requests_.emplace(key, new eckit::mpi::Request);
        //     return collect_map_requests_.at(key);
        // }

    private:
        int this_model_id_;

        std::vector<int> models_;
        std::unordered_map<int, std::string> model2grid_;
        std::unordered_map<int, std::vector<int>> model2ranks_;
        std::unordered_map<std::string, parallel::Collect*> collect_map_;
        // std::unordered_map<std::string, eckit::mpi::Request*> collect_map_requests_;
        // std::unordered_map<int, std::vector<int, ParInter*>> m2m_interpolators_;
    };


    static ModelCoupler coupler_;

    void finalise_coupler();
    void set_default_comm_to_local(int model_id);
    Matrix read_interpolation_matrix(std::string matrix_name);


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
        // parallel::Collect collect;
        if (is_model_2) {
            if (mpi::comm().rank() == 0) {
                gmatrix = read_interpolation_matrix("mat");
            }
            grid2 = atlas::Grid(coupler_.grid_name(model_2));
            grid::Distribution dist_m2(grid2);
            functionspace::StructuredColumns fspace_m2(grid2);

            std::vector<Index> l_rows;
            std::vector<Index> l_cols;
            std::vector<Index> l_gcols;
            std::vector<Value> l_vals;
            ParInter::extract(fspace_m2, dist_m2, gmatrix, l_rows, l_cols, l_gcols, l_vals);

            grid1 = atlas::Grid(coupler_.grid_name(model_1));
            grid::Distribution dist_m1(grid1, grid1.partitioner() | util::Config("nb_partitions", coupler_.global_ranks(model_1).size()));

            std::size_t collect_size_ = l_gcols.size();
            std::vector<gidx_t> collect_gidx(collect_size_);
            std::vector<idx_t>  collect_ridx(collect_size_);
            std::vector<int> collect_partition;
            idx_t ridx_base = 0;

            for (size_t i=0; i < collect_gidx.size(); ++i) {
                collect_gidx[i] = l_gcols[i] + 1;
            }
            collect_partition.resize(collect_size_);
            util::locate_partition(dist_m1, collect_gidx.size(), collect_gidx.data(), collect_partition.data());
            for (int i = 0; i < collect_partition.size(); ++i) {
                auto part = collect_partition[i]; // local rank
                collect_partition[i] = coupler_.global_ranks(model_1)[part]; // "world"-rank
            }
            util::locate_remote_index(mpi::comm("world"), 0, nullptr, nullptr,
                                      collect_size_, collect_gidx.data(), collect_partition.data(), collect_ridx.data(), ridx_base);
            if (mpi::comm().rank() == 0) {
                std::cout << " collect_part " << collect_partition << std::endl;
                std::cout << " collect_ridx " << collect_ridx << std::endl;
            }
            collect->setup("world", collect_size_, collect_partition.data(), collect_ridx.data(), ridx_base);
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


    void setup_coupler(int model_id, Grid grid) {
        set_default_comm_to_local(model_id);
        coupler_.register_model(model_id, grid);
        coupler_.print_models();

        // HACK: we assume only m1 sends field 'f1' to m2
        setup_oneway_remap(10, 21);
    }


    void put_field(Field f, int model_2, int tstep = 0) {
        auto& collect = *(coupler_.collect_map(coupler_.this_model(), model_2));
        collect.send<double, 1>(f);
    }


    void get_field(Field f, int model_1, int tstep = 0) {
        auto& collect = *(coupler_.collect_map(model_1, coupler_.this_model()));
        collect.recv<double, 1>(f);
    }

} // namespace atlas::mct
