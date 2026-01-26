/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
#include "atlas/functionspace/StructuredColumns.h"

#include <functional>
#include <iomanip>
#include <new>  // for bad_alloc exception
#include <sstream>
#include <string>
#include <unordered_map>

#include "eckit/log/Bytes.h"

#include "atlas/array/MakeView.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/library/Library.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/Statistics.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/vector.h"

namespace atlas {
namespace functionspace {
namespace detail {

void StructuredColumns::create_remote_index() const {
    field_remote_index_ = Field("remote_idx", array::make_datatype<idx_t>(), array::make_shape(size_halo_));
    auto remote_idx     = array::make_indexview<idx_t, 1>(field_remote_index_);
    atlas_omp_parallel_for(idx_t n = 0; n < size_owned_; ++n) { remote_idx(n) = n; }


    ATLAS_TRACE_SCOPE("Parallelisation ...") {
        auto build_partition_graph = [this]() -> std::unique_ptr<Mesh::PartitionGraph> {
            const eckit::mpi::Comm& comm = mpi::comm(mpi_comm());
            const int mpi_size           = int(comm.size());
            const int mpi_rank           = int(comm.rank());

            auto part = array::make_view<int, 1>(this->partition());

            std::vector<int> others_set(mpi_size, 0);
            others_set[mpi_rank] = 1;
            for (idx_t i = size_owned_; i < size_halo_; ++i) {
                others_set[part(i)] = 1;  // present
            }
            std::vector<int> others;
            others.reserve(mpi_size);
            for (idx_t p = 0; p < mpi_size; ++p) {
                if (others_set[p]) {
                    others.emplace_back(p);
                }
            }

            eckit::mpi::Buffer<int> recv_others(mpi_size);

            ATLAS_TRACE_MPI(ALLGATHER) { comm.allGatherv(others.begin(), others.end(), recv_others); }

            std::vector<idx_t> counts(recv_others.counts.begin(), recv_others.counts.end());
            std::vector<idx_t> displs(recv_others.displs.begin(), recv_others.displs.end());
            std::vector<idx_t> values(recv_others.buffer.begin(), recv_others.buffer.end());
            return std::unique_ptr<Mesh::PartitionGraph>(
                new Mesh::PartitionGraph(values.data(), mpi_size, displs.data(), counts.data()));
        };

        std::unique_ptr<Mesh::PartitionGraph> graph_ptr;
        ATLAS_TRACE_SCOPE("Building partition graph...") { graph_ptr = build_partition_graph(); }
        const Mesh::PartitionGraph& graph = *graph_ptr;

        ATLAS_TRACE_SCOPE("Setup remote_index fields...") {
            auto p = array::make_view<int, 1>(partition());
            auto g = array::make_view<gidx_t, 1>(global_index());

            const eckit::mpi::Comm& comm = mpi::comm(mpi_comm());
            const int mpi_rank           = int(comm.rank());

            auto neighbours           = graph.nearestNeighbours(mpi_rank);
            const idx_t nb_neighbours = static_cast<idx_t>(neighbours.size());
            std::unordered_map<int, idx_t> part_to_neighbour;
            ATLAS_TRACE_SCOPE("part_to_neighbour") {
                part_to_neighbour.reserve(nb_neighbours);
                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    part_to_neighbour[neighbours[j]] = j;
                }
            }
            std::vector<idx_t> halo_per_neighbour(neighbours.size(), 0);
            ATLAS_TRACE_SCOPE("set halo_per_neighbour")
            for (idx_t i = size_owned_; i < size_halo_; ++i) {
                halo_per_neighbour[part_to_neighbour[p(i)]]++;
            }

            std::vector<std::vector<gidx_t>> g_per_neighbour(neighbours.size());
            ATLAS_TRACE_SCOPE("assemble g_per_neighbour") {
                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    g_per_neighbour[j].reserve(halo_per_neighbour[j]);
                }
                for (idx_t j = size_owned_; j < size_halo_; ++j) {
                    g_per_neighbour[part_to_neighbour[p(j)]].emplace_back(g(j));
                }
            }

            std::vector<eckit::mpi::Request> send_requests(neighbours.size());
            std::vector<eckit::mpi::Request> recv_requests(neighbours.size());

            std::vector<idx_t> recv_size(neighbours.size());
            std::vector<idx_t> send_size(neighbours.size());

            int tag = 0;
            ATLAS_TRACE_SCOPE("send-receive g_per_neighbour size") {
                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    send_size[j]     = static_cast<idx_t>(g_per_neighbour[j].size());
                    send_requests[j] = comm.iSend(send_size[j], neighbours[j], tag);
                    recv_requests[j] = comm.iReceive(recv_size[j], neighbours[j], tag);
                }

                ATLAS_TRACE_MPI(WAIT) {
                    for (idx_t j = 0; j < nb_neighbours; ++j) {
                        comm.wait(send_requests[j]);
                    }

                    for (idx_t j = 0; j < nb_neighbours; ++j) {
                        comm.wait(recv_requests[j]);
                    }
                }
            }

            std::vector<std::vector<gidx_t>> recv_g_per_neighbour(neighbours.size());
            ATLAS_TRACE_SCOPE("send-receive g_per_neighbour")
            for (idx_t j = 0; j < nb_neighbours; ++j) {
                recv_g_per_neighbour[j].resize(recv_size[j]);

                send_requests[j] = comm.iSend(g_per_neighbour[j].data(), g_per_neighbour[j].size(), neighbours[j], tag);
                recv_requests[j] =
                    comm.iReceive(recv_g_per_neighbour[j].data(), recv_g_per_neighbour[j].size(), neighbours[j], tag);
            }

            std::vector<std::vector<idx_t>> send_r_per_neighbour(neighbours.size());

            // Assemble "send_r_per_neighbour"
            // Depending if we can afford memory for a globally sized array, we can have a
            // much faster version of g_to_r map using std::vector.
            // TODO: using c++14 we can create a polymorphic lambda to avoid duplicated
            // code in the two branches of following if.

            size_t max_glb_idx = grid().size();
            atlas::vector<idx_t> g_to_r_vector;
            bool use_unordered_map_fallback = false;
            try {
                g_to_r_vector.resize(max_glb_idx + 1);
            }
            catch (std::bad_alloc& e) {
                if (comm.size() > 1) {
                    Log::warning() << "Could not allocate " << eckit::Bytes((max_glb_idx + 1) * sizeof(idx_t)) << Here()
                                   << "\n"
                                   << "Using slower unordered_map fallback to map global to remote indices"
                                   << std::endl;
                    use_unordered_map_fallback = true;
                }
                else {
                    throw_Exception(
                        "Could not allocate " + std::string(eckit::Bytes((max_glb_idx + 1) * sizeof(idx_t))), Here());
                }
            }
            if (not use_unordered_map_fallback) {
                auto& g_to_r = g_to_r_vector;
                ATLAS_TRACE_SCOPE("g_to_r (using vector)") {
                    atlas_omp_parallel_for(idx_t j = 0; j < size_owned_; ++j) {
                        // parallel omp possible for ``` g_to_r[g(j)] = j ``` as we only loop over size_owned,
                        // where global_index is such that race-conditions cannot occur
#if ATLAS_ARRAYVIEW_BOUNDS_CHECKING
                        if (g(j) >= max_glb_idx + 1) {
                            ATLAS_DEBUG_VAR(g(j));
                            throw_OutOfRange("g_to_r", g(j), max_glb_idx + 1, Here());
                        }
#endif
                        g_to_r[g(j)] = j;
                    }
                }
                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    send_r_per_neighbour[j].reserve(recv_size[j]);

                    comm.wait(recv_requests[j]);  // wait for recv_g_per_neighbour[j]
                    for (gidx_t gidx : recv_g_per_neighbour[j]) {
                        send_r_per_neighbour[j].emplace_back(g_to_r[gidx]);
                    }
                }
            }
            else {
                std::unordered_map<gidx_t, idx_t> g_to_r;
                g_to_r.reserve(size_owned_);

                ATLAS_TRACE_SCOPE("g_to_r (using unordered set)") {
                    for (idx_t j = 0; j < size_owned_; ++j) {
                        g_to_r[g(j)] = j;
                    }
                }
                ATLAS_TRACE_SCOPE("assemble r_per_neighbour")
                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    send_r_per_neighbour[j].reserve(recv_size[j]);

                    comm.wait(recv_requests[j]);  // wait for recv_g_per_neighbour[j]
                    for (gidx_t gidx : recv_g_per_neighbour[j]) {
                        send_r_per_neighbour[j].emplace_back(g_to_r[gidx]);
                    }
                }
            }

            std::vector<std::vector<idx_t>> r_per_neighbour(neighbours.size());

            ATLAS_TRACE_SCOPE("send-receive r_per_neighbour") {
                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    r_per_neighbour[j].resize(halo_per_neighbour[j]);
                }

                ATLAS_TRACE_MPI(ALLTOALL) {
                    for (idx_t j = 0; j < nb_neighbours; ++j) {
                        comm.wait(send_requests[j]);
                        send_requests[j] = comm.iSend(send_r_per_neighbour[j].data(), send_r_per_neighbour[j].size(),
                                                      neighbours[j], tag);
                        recv_requests[j] =
                            comm.iReceive(r_per_neighbour[j].data(), r_per_neighbour[j].size(), neighbours[j], tag);
                    }
                }
                ATLAS_TRACE_MPI(WAIT) {
                    for (idx_t j = 0; j < nb_neighbours; ++j) {
                        comm.wait(recv_requests[j]);
                    }
                }
            }


            std::vector<idx_t> counters(neighbours.size(), 0);
            ATLAS_TRACE_SCOPE("set remote_idx")
            for (idx_t j = size_owned_; j < size_halo_; ++j) {
                idx_t neighbour = part_to_neighbour[p(j)];
                remote_idx(j)   = r_per_neighbour[neighbour][counters[neighbour]++];
            }

            ATLAS_TRACE_MPI(WAIT) {
                for (idx_t j = 0; j < nb_neighbours; ++j) {
                    comm.wait(send_requests[j]);
                }
            }
        }
    }
}

}  // namespace detail

// ----------------------------------------------------------------------------

}  // namespace functionspace
}  // namespace atlas
