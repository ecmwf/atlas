/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <iostream>

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildPeriodicBoundaries.h"
#include "atlas/parallel/mpi/Statistics.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/PeriodicTransform.h"

using Topology = atlas::mesh::Nodes::Topology;
using atlas::util::LonLatMicroDeg;
using atlas::util::PeriodicTransform;

namespace atlas {
namespace mesh {
namespace actions {

using uid_t = gidx_t;

void build_periodic_boundaries(Mesh& mesh) {
    ATLAS_TRACE();
    bool periodic = false;
    mesh.metadata().get("periodic", periodic);
    if (!periodic) {
        mpi::Scope mpi_scope(mesh.mpi_comm());

        auto mpi_size = mpi::size();
        auto mypart   = mpi::rank();

        mesh::Nodes& nodes = mesh.nodes();

        auto flags = array::make_view<int, 1>(nodes.flags());
        auto ridx  = array::make_indexview<idx_t, 1>(nodes.remote_index());
        auto part  = array::make_view<int, 1>(nodes.partition());
        auto ghost = array::make_view<int, 1>(nodes.ghost());

        int nb_nodes = nodes.size();

        auto xy = array::make_view<double, 2>(nodes.xy());

        // Identify my master and slave nodes on own partition
        // master nodes are at x=0,  slave nodes are at x=2pi
        std::map<uid_t, int> master_lookup;
        std::map<uid_t, int> slave_lookup;
        std::vector<int> master_nodes;
        master_nodes.reserve(3 * nb_nodes);
        std::vector<int> slave_nodes;
        slave_nodes.reserve(3 * nb_nodes);

        for (idx_t jnode = 0; jnode < nodes.size(); ++jnode) {
            if (Topology::check_all(flags(jnode), Topology::BC | Topology::WEST)) {
                Topology::set(flags(jnode), Topology::PERIODIC);
                if (part(jnode) == mypart) {
                    LonLatMicroDeg ll(xy(jnode, XX), xy(jnode, YY));
                    master_lookup[ll.unique()] = jnode;
                    master_nodes.push_back(ll.lon());
                    master_nodes.push_back(ll.lat());
                    master_nodes.push_back(jnode);
                }
            }
            else if (Topology::check(flags(jnode), Topology::BC | Topology::EAST)) {
                Topology::set(flags(jnode), Topology::PERIODIC);
                Topology::set(flags(jnode), Topology::GHOST);
                ghost(jnode) = 1;
                LonLatMicroDeg ll(xy(jnode, XX), xy(jnode, YY));
                slave_lookup[ll.unique()] = jnode;
                slave_nodes.push_back(ll.lon());
                slave_nodes.push_back(ll.lat());
                slave_nodes.push_back(jnode);
                ridx(jnode) = -1;
            }
        }

        std::vector<std::vector<int>> found_master(mpi_size);
        std::vector<std::vector<int>> send_slave_idx(mpi_size);

        // Find masters on other tasks to send to me
        {
            int sendcnt = slave_nodes.size();
            std::vector<int> recvcounts(mpi_size);

            ATLAS_TRACE_MPI(ALLGATHER) { mpi::comm().allGather(sendcnt, recvcounts.begin(), recvcounts.end()); }

            std::vector<int> recvdispls(mpi_size);
            recvdispls[0] = 0;
            int recvcnt   = recvcounts[0];
            for (idx_t jproc = 1; jproc < mpi_size; ++jproc) {
                recvdispls[jproc] = recvdispls[jproc - 1] + recvcounts[jproc - 1];
                recvcnt += recvcounts[jproc];
            }
            std::vector<int> recvbuf(recvcnt);

            ATLAS_TRACE_MPI(ALLGATHER) {
                mpi::comm().allGatherv(slave_nodes.begin(), slave_nodes.end(), recvbuf.begin(), recvcounts.data(),
                                       recvdispls.data());
            }

            PeriodicTransform transform;
            for (idx_t jproc = 0; jproc < mpi_size; ++jproc) {
                found_master.reserve(master_nodes.size());
                send_slave_idx.reserve(master_nodes.size());
                array::LocalView<int, 2> recv_slave(recvbuf.data() + recvdispls[jproc],
                                                    array::make_shape(recvcounts[jproc] / 3, 3));
                for (idx_t jnode = 0; jnode < recv_slave.shape(0); ++jnode) {
                    LonLatMicroDeg slave(recv_slave(jnode, LON), recv_slave(jnode, LAT));
                    transform(slave, -1);
                    uid_t slave_uid = slave.unique();
                    if (master_lookup.count(slave_uid)) {
                        int master_idx = master_lookup[slave_uid];
                        int slave_idx  = recv_slave(jnode, 2);
                        found_master[jproc].push_back(master_idx);
                        send_slave_idx[jproc].push_back(slave_idx);
                    }
                }
            }
        }

        // Fill in data to communicate
        std::vector<std::vector<int>> recv_slave_idx(mpi_size);
        std::vector<std::vector<int>> send_master_part(mpi_size);
        std::vector<std::vector<int>> recv_master_part(mpi_size);
        std::vector<std::vector<int>> send_master_ridx(mpi_size);
        std::vector<std::vector<int>> recv_master_ridx(mpi_size);

        //  std::vector< std::vector<int> > send_slave_part( mpi_size );
        //  std::vector< std::vector<int> > recv_slave_part( mpi_size );
        //  std::vector< std::vector<int> > send_slave_ridx( mpi_size );
        //  std::vector< std::vector<int> > recv_slave_ridx( mpi_size );

        {
            for (idx_t jproc = 0; jproc < mpi_size; ++jproc) {
                idx_t nb_found_master = static_cast<idx_t>(found_master[jproc].size());
                send_master_part[jproc].resize(nb_found_master);
                send_master_ridx[jproc].resize(nb_found_master);
                for (idx_t jnode = 0; jnode < nb_found_master; ++jnode) {
                    int loc_idx                    = found_master[jproc][jnode];
                    send_master_part[jproc][jnode] = part(loc_idx);
                    send_master_ridx[jproc][jnode] = loc_idx;
                }

                //      int nb_found_slaves = found_slave[jproc].size();
                //      send_slave_glb_idx[jproc].resize(nb_found_slaves);
                //      send_slave_part   [jproc].resize(nb_found_slaves);
                //      send_slave_ridx   [jproc].resize(nb_found_slaves);
                //      for( int jnode=0; jnode<nb_found_slaves; ++jnode )
                //      {
                //        int loc_idx = found_slave[jproc][jnode];
                //        send_slave_glb_idx[jproc][jnode] = glb_idx( loc_idx );
                //        send_slave_part   [jproc][jnode] = part   ( loc_idx );
                //        send_slave_ridx   [jproc][jnode] = loc_idx;
                //      }
            }
        }

        // Communicate
        ATLAS_TRACE_MPI(ALLTOALL) {
            mpi::comm().allToAll(send_slave_idx, recv_slave_idx);
            mpi::comm().allToAll(send_master_part, recv_master_part);
            mpi::comm().allToAll(send_master_ridx, recv_master_ridx);
            //  mpi::comm().allToAll( send_slave_part,     recv_slave_part    );
            //  mpi::comm().allToAll( send_slave_loc,      recv_slave_ridx    );
        }

        // Fill in periodic
        // unused // int nb_recv_master = 0;
        for (idx_t jproc = 0; jproc < mpi_size; ++jproc) {
            idx_t nb_recv = static_cast<idx_t>(recv_slave_idx[jproc].size());
            for (idx_t jnode = 0; jnode < nb_recv; ++jnode) {
                idx_t slave_idx = recv_slave_idx[jproc][jnode];
                part(slave_idx) = recv_master_part[jproc][jnode];
                ridx(slave_idx) = recv_master_ridx[jproc][jnode];
            }
        }
    }
    mesh.metadata().set("periodic", true);
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_periodic_boundaries(Mesh::Implementation* mesh) {
    ATLAS_ASSERT(mesh != nullptr, "Cannot access uninitialised atlas_Mesh");
    Mesh m(mesh);
    build_periodic_boundaries(m);
}
// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
