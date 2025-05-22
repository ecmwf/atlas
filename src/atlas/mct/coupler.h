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
        struct {
            // std::vector<eckit::mpi::Request> recv_req;
            size_t src_size;
            size_t src_field_shape;
            parallel::Collect* collect;
        } field_exchange;

    private:
        int this_model_id_;

        std::vector<int> models_;
        std::unordered_map<int, std::string> model2grid_;
        std::unordered_map<int, std::vector<int>> model2ranks_;
        std::unordered_map<std::string, parallel::Collect*> collect_map_;
        // std::unordered_map<int, std::vector<int, ParInter*>> m2m_fields_;
    };

    static ModelCoupler coupler_;

} // namespace atlas::mct
