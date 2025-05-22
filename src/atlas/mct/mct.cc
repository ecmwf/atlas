/*
 * (C) Copyright 2025 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mct/mct.h"

using Matrix      = atlas::linalg::SparseMatrixStorage;

namespace atlas::mct {

    void ModelCoupler::register_model(int model_id, Grid grid) {
        std::vector<int> models;
        auto& gcomm = mpi::comm("world");
        auto gcommsize = gcomm.size();

        // exchange models
        models.resize(gcommsize);
        models[gcomm.rank()] = model_id;
        gcomm.allGather(model_id, models.begin(), models.end());

        // exchange other procs' global ranks
        std::vector<int> commsizes(gcommsize);
        int my_size = mpi::comm().size();
        gcomm.allGather(my_size, commsizes.begin(), commsizes.end());

        std::vector<int> my_global_ranks(my_size);
        int my_global_rank = gcomm.rank();
        mpi::comm().allGather(my_global_rank, my_global_ranks.begin(), my_global_ranks.end());

        eckit::mpi::Buffer<int> buf_ranks(gcommsize);
        gcomm.allGatherv(my_global_ranks.begin(), my_global_ranks.end(), buf_ranks);

        // exchange grid names
        std::vector<char> gn(grid.name().size());
        for (int i = 0; i < grid.name().size(); i++) {
            gn[i] = grid.name()[i];
        }
        eckit::mpi::Buffer<char> buf(gcommsize);
        gcomm.allGatherv(gn.begin(), gn.end(), buf);

        for (int i = 0; i < gcommsize; i++) {
            if (std::find(models_.begin(), models_.end(), models[i]) == models_.end()) {
                models_.emplace_back(models[i]);
                model2grid_[models[i]] = std::string(&(buf.buffer[buf.displs[i]]), buf.counts[i]);
                model2ranks_[models[i]].reserve(buf_ranks.counts[i]);
                std::copy_n(&(buf_ranks.buffer[buf_ranks.displs[i]]), buf_ranks.counts[i], 
                    std::back_inserter(model2ranks_[models[i]]));
            }
        }

        // find my model id
        for (int i = 0; i < models_.size(); ++i) {
            if (mpi::comm().name() == comm_name(models_[i])) {
                this_model_id_ = models_[i];
                break;
            }
        }
    }


    void set_default_comm_to_local(int model_id)
    {
        eckit::mpi::comm().split(model_id, coupler_.comm_name(model_id));
        eckit::mpi::setCommDefault(coupler_.comm_name(model_id));
    }


    Matrix read_interpolation_matrix(std::string matrix_name) {
        eckit::linalg::SparseMatrix eckit_matrix;
        Log::info() << "reading matrix '" << matrix_name << "'." <<std::endl;
        eckit_matrix.load(matrix_name + ".eckit");
        return atlas::linalg::make_sparse_matrix_storage(std::move(eckit_matrix));
    }


    void ModelCoupler::print_models() const {
        std::vector<int> lead_ranks(models_.size());
        for (int i = 0; i < models_.size(); i++) {
            lead_ranks[i] = model2ranks_.at(models_[i])[0];
        }
        for (int i = 0; i < models_.size(); i++) {
            mpi::comm("world").barrier();
            if (mpi::comm("world").rank() == lead_ranks[i]) {
                std::cout << "Model " << models_[i] << ", grid " << model2grid_.at(models_[i])
                    << ", global ranks are : " << model2ranks_.at(models_[i]) << std::endl;
            }
        }
    }


    void ModelCoupler::finalise() {
        models_.clear();
        model2grid_.clear();
        model2ranks_.clear();
    }


    void finalise_coupler() {
        coupler_.finalise();
    }

} // namespace atlas::mct
