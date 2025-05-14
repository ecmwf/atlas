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
        void register_model(int model_id, Grid grid) {
            models_.emplace_back(model_id);
            grids_.emplace_back(grid);
        }
        static std::string comm_name(int model_id) { return std::to_string(model_id); }

        Grid get_grid(int model_id) {
            auto it = std::find(models_.begin(), models_.end(), model_id);
            return grids_[std::distance(models_.begin(), it)];
        }

    private:
        std::vector<int> models_;
        std::vector<Grid> grids_;
    };

    static ModelCoupler coupler_;


    Matrix read_interpolation_matrix(std::string matrix_name) {
        eckit::linalg::SparseMatrix eckit_matrix;
        Log::info() << "reading matrix '" << matrix_name << "'." <<std::endl;
        eckit_matrix.load(matrix_name + ".eckit");
        return atlas::linalg::make_sparse_matrix_storage(std::move(eckit_matrix));
    }


    void set_default_comm_to_local(int model_id) {
        eckit::mpi::comm().split(model_id, ModelCoupler::comm_name(model_id));
        eckit::mpi::setCommDefault(ModelCoupler::comm_name(model_id));
    }

    void setup_coupler(int model_id, Grid grid) {
        coupler_.register_model(model_id, grid);
        set_default_comm_to_local(model_id);
        auto matrix = read_interpolation_matrix("mat");

        // Hack: we assume m1 sends field 'f1' to m2
        functionspace::StructuredColumns fspace_m1(grid);
        functionspace::StructuredColumns fspace_m2(grid);
        interpolation::distribute_global_matrix(fspace_m1, fspace_m2, matrix);
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
