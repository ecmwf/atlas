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
#include "atlas/parallel/mpi/mpi.h"

using Matrix      = atlas::linalg::SparseMatrixStorage;

namespace atlas::mct {

    class ModelCoupler {
    public:
        ModelCoupler() {}

        void register_model(int model_id, Grid grid);

        void print_models() const;

        std::vector<int> models() const { return models_; }

        static std::string comm_name(int model_id) { return std::to_string(model_id); }

        Grid get_grid(int model_id) { return model2grid_[model_id]; }

    private:
        std::vector<int> models_;
        std::unordered_map<int, std::string> model2grid_;
        std::unordered_map<int, std::vector<int>> model2ranks_;
    };

    static ModelCoupler coupler_;


    Matrix read_interpolation_matrix(std::string matrix_name) {
        eckit::linalg::SparseMatrix eckit_matrix;
        Log::info() << "reading matrix '" << matrix_name << "'." <<std::endl;
        eckit_matrix.load(matrix_name + ".eckit");
        return atlas::linalg::make_sparse_matrix_storage(std::move(eckit_matrix));
    }


    void set_default_comm_to_local(int model_id)
    {
        eckit::mpi::comm().split(model_id, ModelCoupler::comm_name(model_id));
        eckit::mpi::setCommDefault(ModelCoupler::comm_name(model_id));
    }


    void setup_oneway_remap(int model_1, int model_2) {
        bool is_model_1 = (mpi::comm().name() == ModelCoupler::comm_name(model_1));
        bool is_model_2 = (mpi::comm().name() == ModelCoupler::comm_name(model_2));

        // everything is set up on model_2 for remapping
        // global matrix is read in on model_2
        // field 1 is send to model_2
        if (is_model_2) {
            auto matrix = read_interpolation_matrix("mat");
            functionspace::StructuredColumns fspace_m1(coupler_.get_grid(model_1));
            functionspace::StructuredColumns fspace_m2(coupler_.get_grid(model_2));
            auto lmatrix = interpolation::distribute_global_matrix(fspace_m1, fspace_m2, matrix);
        }
    }


    void setup_coupler(int model_id, Grid grid) {
        set_default_comm_to_local(model_id);
        coupler_.register_model(model_id, grid);

        if (mpi::comm("world").rank() == 5) coupler_.print_models();

        // HACK: we assume only m1 sends field 'f1' to m2
        setup_oneway_remap(10, 21);
    }


    void finalise_coupler() {
        // nothing yet
    }

    void put_field(Field f) {
        // nothing yet
    }

    void get_field(Field f) {
        // nothing yet
    }

} // namespace atlas::mct
