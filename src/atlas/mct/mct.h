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

#include "atlas/grid.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/interpolation/AssembleGlobalMatrix.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/mct/ParInter.h"
#include "atlas/parallel/mpi/mpi.h"

using Matrix      = atlas::linalg::SparseMatrixStorage;

namespace atlas::mct {

    class ModelCoupler {
    public:
        ModelCoupler() {}

        void register_model(int model_id, Grid grid);
        void print_models() const;

        const std::vector<int>& models() const { return models_; }
        const std::string comm_name(int model_id) const { return std::to_string(model_id); }

        const std::string& grid_name(int model_id) const { return model2grid_.at(model_id); }
        const std::vector<int>& global_ranks(int model_id) const { return model2ranks_.at(model_id); }

        void finalise() {
            models_.clear();
            model2grid_.clear();
            model2ranks_.clear();
        }

    private:
        std::vector<int> models_;
        std::unordered_map<int, std::string> model2grid_;
        std::unordered_map<int, std::vector<int>> model2ranks_;
    };


    static ModelCoupler coupler_;

    void finalise_coupler();
    void set_default_comm_to_local(int model_id);
    Matrix read_interpolation_matrix(std::string matrix_name);


    void setup_oneway_remap(int model_1, int model_2) {
        bool is_model_1 = (mpi::comm().name() == coupler_.comm_name(model_1));
        bool is_model_2 = (mpi::comm().name() == coupler_.comm_name(model_2));

        // everything is set up on model_2 for remapping
        // global matrix is read in on model_2
        // field 1 is send to model_2

        atlas::Grid grid1;
        atlas::Grid grid2; // later no need to recreate grid2, as we are on the tasks which own grid2
        FunctionSpace fspace_m1;
        FunctionSpace fspace_m2;

        Matrix gmatrix;
        if (is_model_2 && mpi::comm().rank() == 0) {
            gmatrix = read_interpolation_matrix("mat");
        }

        // prepare the communication pattern for sending field data from model_1 to model_2
        if (is_model_1 || is_model_2) {
            grid1 = atlas::Grid(coupler_.grid_name(model_1));
            grid2 = atlas::Grid(coupler_.grid_name(model_2));
            fspace_m1 = functionspace::StructuredColumns(grid1);
            fspace_m2 = functionspace::StructuredColumns(grid2);

            auto partitioner_m1_config = grid1.partitioner();
            partitioner_m1_config.set("partitions", coupler_.global_ranks(model_1).size());
            grid::Distribution dist_m1( grid1, grid::Partitioner(partitioner_m1_config) );
            // map dist_m1 to ranks in world communicator

            auto partitioner_m2_config = grid1.partitioner();
            partitioner_m2_config.set("partitions", coupler_.global_ranks(model_2).size());
            grid::Distribution dist_m2( grid1, grid::Partitioner(partitioner_m2_config) );

            // ParInter(gmatrix, fspace_m1, fspace_m2, dist_m1, dist_m2);
        }

        Matrix lmatrix;
        if (is_model_2) {
            lmatrix = interpolation::distribute_global_matrix(fspace_m1, fspace_m2, gmatrix);
        }
    }


    void setup_coupler(int model_id, Grid grid) {
        set_default_comm_to_local(model_id);
        coupler_.register_model(model_id, grid);

        coupler_.print_models();

        // HACK: we assume only m1 sends field 'f1' to m2
        setup_oneway_remap(10, 21);
    }


    void put_field(Field f) {
        // nothing yet
    }


    void get_field(Field f) {
        // nothing yet
    }

} // namespace atlas::mct
