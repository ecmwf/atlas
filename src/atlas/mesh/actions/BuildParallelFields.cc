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
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/Buffer.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/PeriodicTransform.h"
#include "atlas/util/Unique.h"

#define EDGE(jedge)                                                                               \
    "Edge(" << node_gidx(edge_nodes(jedge, 0)) << "[p" << node_part(edge_nodes(jedge, 0)) << "] " \
            << node_gidx(edge_nodes(jedge, 1)) << "[p" << node_part(edge_nodes(jedge, 1)) << "])"

//#define DEBUGGING_PARFIELDS
#ifdef DEBUGGING_PARFIELDS
#define own1 2419089
#define own2 2423185
#define OWNED_EDGE(jedge)                                                                    \
    ((node_gidx(edge_nodes(jedge, 0)) == own1 && node_gidx(edge_nodes(jedge, 1)) == own2) || \
     (node_gidx(edge_nodes(jedge, 0)) == own2 && node_gidx(edge_nodes(jedge, 1)) == own1))
#define per1 -1
#define per2 -1
#define PERIODIC_EDGE(jedge)                                                                 \
    ((node_gidx(edge_nodes(jedge, 0)) == per1 && node_gidx(edge_nodes(jedge, 1)) == per2) || \
     (node_gidx(edge_nodes(jedge, 0)) == per2 && node_gidx(edge_nodes(jedge, 1)) == per1))
#define find1 39  // 1697
#define find2 41  // 1698
#define FIND_EDGE(jedge)                                                                       \
    ((node_gidx(edge_nodes(jedge, 0)) == find1 && node_gidx(edge_nodes(jedge, 1)) == find2) || \
     (node_gidx(edge_nodes(jedge, 0)) == find2 && node_gidx(edge_nodes(jedge, 1)) == find1))
#define find_gidx 689849552510167040
#define FIND_GIDX(UID) ((UID) == find_gidx)
#endif

using Topology = atlas::mesh::Nodes::Topology;
using atlas::util::PeriodicTransform;
using atlas::util::UniqueLonLat;

namespace atlas {
namespace mesh {
namespace actions {

Field& build_nodes_partition(mesh::Nodes& nodes);
Field& build_nodes_remote_idx(mesh::Nodes& nodes);
Field& build_nodes_global_idx(mesh::Nodes& nodes);
Field& build_edges_partition(Mesh& mesh);
Field& build_edges_remote_idx(Mesh& mesh);
Field& build_edges_global_idx(Mesh& mesh);

//----------------------------------------------------------------------------------------------------------------------

using uid_t = gidx_t;

namespace {

struct Node {
    Node(gidx_t gid, idx_t idx) {
        g = gid;
        i = idx;
    }
    gidx_t g;
    idx_t i;
    bool operator<(const Node& other) const { return (g < other.g); }
};

}  // namespace

//----------------------------------------------------------------------------------------------------------------------

void build_parallel_fields(Mesh& mesh) {
    ATLAS_TRACE();
    build_nodes_parallel_fields(mesh);
}

//----------------------------------------------------------------------------------------------------------------------

void build_nodes_parallel_fields(Mesh& mesh) {
    mpi::Scope mpi_scope(mesh.mpi_comm());
    build_nodes_parallel_fields(mesh.nodes());
}

void build_nodes_parallel_fields(mesh::Nodes& nodes) {
    ATLAS_TRACE();
    bool parallel = false;
    nodes.metadata().get("parallel", parallel);
    if (!parallel) {
        build_nodes_partition(nodes);
        build_nodes_remote_idx(nodes);
        build_nodes_global_idx(nodes);
    }
    nodes.metadata().set("parallel", true);
}

//----------------------------------------------------------------------------------------------------------------------

void build_edges_parallel_fields(Mesh& mesh) {
    ATLAS_TRACE();
    mpi::Scope mpi_scope(mesh.mpi_comm());
    build_edges_partition(mesh);
    build_edges_remote_idx(mesh);
    /*
 * We turn following off. It is expensive and we don't really care about a nice
 * contiguous
 * ordering.
 */
    // build_edges_global_idx( mesh );
}

//----------------------------------------------------------------------------------------------------------------------

Field& build_nodes_global_idx(mesh::Nodes& nodes) {
    ATLAS_TRACE();

    array::ArrayView<gidx_t, 1> glb_idx = array::make_view<gidx_t, 1>(nodes.global_index());

    UniqueLonLat compute_uid(nodes);

    for (idx_t jnode = 0; jnode < glb_idx.shape(0); ++jnode) {
        if (glb_idx(jnode) <= 0) {
            glb_idx(jnode) = compute_uid(jnode);
        }
    }
    return nodes.global_index();
}

void renumber_nodes_glb_idx(Mesh& mesh) {
    mpi::Scope mpi_scope(mesh.mpi_comm());
    renumber_nodes_glb_idx(mesh.nodes());
}

void renumber_nodes_glb_idx(mesh::Nodes& nodes) {
    bool human_readable(false);
    nodes.global_index().metadata().get("human_readable", human_readable);
    if (human_readable) {
        /* nothing to be done */
        return;
    }

    ATLAS_TRACE();

    // TODO: ATLAS-14: fix renumbering of EAST periodic boundary points
    // --> Those specific periodic points at the EAST boundary are not checked for
    // uid,
    //     and could receive different gidx for different tasks

    UniqueLonLat compute_uid(nodes);

    // unused // int mypart = mpi::rank();
    int nparts = mpi::size();
    idx_t root = 0;

    array::ArrayView<gidx_t, 1> glb_idx = array::make_view<gidx_t, 1>(nodes.global_index());

    /*
 * Sorting following gidx will define global order of
 * gathered fields. Special care needs to be taken for
 * pole edges, as their centroid might coincide with
 * other edges
 */
    int nb_nodes = glb_idx.shape(0);
    for (int jnode = 0; jnode < nb_nodes; ++jnode) {
        if (glb_idx(jnode) <= 0) {
            glb_idx(jnode) = compute_uid(jnode);
        }
    }

    // 1) Gather all global indices, together with location
    array::ArrayT<uid_t> loc_id_arr(nb_nodes);
    array::ArrayView<uid_t, 1> loc_id = array::make_view<uid_t, 1>(loc_id_arr);

    for (int jnode = 0; jnode < nb_nodes; ++jnode) {
        loc_id(jnode) = glb_idx(jnode);
    }

    std::vector<int> recvcounts(mpi::size());
    std::vector<int> recvdispls(mpi::size());

    ATLAS_TRACE_MPI(GATHER) { mpi::comm().gather(nb_nodes, recvcounts, root); }

    recvdispls[0] = 0;
    for (int jpart = 1; jpart < nparts; ++jpart) {  // start at 1
        recvdispls[jpart] = recvcounts[jpart - 1] + recvdispls[jpart - 1];
    }
    int glb_nb_nodes = std::accumulate(recvcounts.begin(), recvcounts.end(), 0);

    array::ArrayT<uid_t> glb_id_arr(glb_nb_nodes);
    array::ArrayView<uid_t, 1> glb_id = array::make_view<uid_t, 1>(glb_id_arr);

    ATLAS_TRACE_MPI(GATHER) {
        mpi::comm().gatherv(loc_id.data(), loc_id.size(), glb_id.data(), recvcounts.data(), recvdispls.data(), root);
    }

    // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
    std::vector<Node> node_sort;
    node_sort.reserve(glb_nb_nodes);
    ATLAS_TRACE_SCOPE("sort global indices") {
        for (idx_t jnode = 0; jnode < glb_id.shape(0); ++jnode) {
            node_sort.emplace_back(glb_id(jnode), jnode);
        }
        std::sort(node_sort.begin(), node_sort.end());
    }

    // Assume edge gid start
    uid_t gid                   = 0;
    const idx_t nb_sorted_nodes = static_cast<idx_t>(node_sort.size());
    for (idx_t jnode = 0; jnode < nb_sorted_nodes; ++jnode) {
        if (jnode == 0) {
            ++gid;
        }
        else if (node_sort[jnode].g != node_sort[jnode - 1].g) {
            ++gid;
        }
        idx_t inode   = node_sort[jnode].i;
        glb_id(inode) = gid;
    }

    // 3) Scatter renumbered back
    ATLAS_TRACE_MPI(SCATTER) {
        mpi::comm().scatterv(glb_id.data(), recvcounts.data(), recvdispls.data(), loc_id.data(), loc_id.size(), root);
    }

    for (int jnode = 0; jnode < nb_nodes; ++jnode) {
        glb_idx(jnode) = loc_id(jnode);
    }
    nodes.global_index().metadata().set("human_readable", true);
}

//----------------------------------------------------------------------------------------------------------------------

Field& build_nodes_remote_idx(mesh::Nodes& nodes) {
    ATLAS_TRACE();
    idx_t mypart = static_cast<idx_t>(mpi::rank());
    idx_t nparts = static_cast<idx_t>(mpi::size());


    std::vector<idx_t> proc(nparts);
    for (idx_t jpart = 0; jpart < nparts; ++jpart) {
        proc[jpart] = jpart;
    }

    auto ridx         = array::make_indexview<idx_t, 1>(nodes.remote_index());
    const auto part   = array::make_view<int, 1>(nodes.partition());
    const auto flags  = array::make_view<int, 1>(nodes.flags());
    const auto lonlat = array::make_view<double, 2>(nodes.lonlat());

    const PeriodicTransform transform_periodic_east(-360.);
    const PeriodicTransform transform_periodic_west(+360.);
    const UniqueLonLat compute_uid_lonlat(nodes);

    auto compute_uid = [&](idx_t jnode) {
        constexpr int PERIODIC = util::Topology::PERIODIC;
        constexpr int EAST     = util::Topology::EAST;
        constexpr int WEST     = util::Topology::WEST;
        const auto flags_view  = util::Bitflags::view(flags(jnode));
        if (flags_view.check(PERIODIC | EAST)) {
            return compute_uid_lonlat(jnode, transform_periodic_east);
        }
        if (flags_view.check(PERIODIC | WEST)) {
            return compute_uid_lonlat(jnode, transform_periodic_west);
        }
        return compute_uid_lonlat(jnode);
    };

    idx_t nb_nodes = nodes.size();

    constexpr idx_t varsize = 2;

    std::vector<std::vector<uid_t>> send_needed(mpi::size());
    std::vector<std::vector<uid_t>> recv_needed(mpi::size());
    int sendcnt = 0;
    std::map<uid_t, int> lookup;
    for (idx_t jnode = 0; jnode < nb_nodes; ++jnode) {
        uid_t uid = compute_uid(jnode);

        if (idx_t(part(jnode)) == mypart) {
            lookup[uid] = jnode;
            ridx(jnode) = jnode;
        }
        else {
            ATLAS_ASSERT(jnode < part.shape(0));
            if (part(jnode) >= static_cast<int>(proc.size())) {
                std::stringstream msg;
                msg << "Assertion [part(" << jnode << ") < proc.size()] failed\n"
                    << "part(" << jnode << ") = " << part(jnode) << "\n"
                    << "proc.size() = " << proc.size();
                throw_AssertionFailed(msg.str(), Here());
            }
            ATLAS_ASSERT(part(jnode) < (idx_t)proc.size());
            ATLAS_ASSERT((size_t)proc[part(jnode)] < send_needed.size());
            send_needed[proc[part(jnode)]].push_back(uid);
            send_needed[proc[part(jnode)]].push_back(jnode);
            sendcnt++;
        }
    }

    ATLAS_TRACE_MPI(ALLTOALL) { mpi::comm().allToAll(send_needed, recv_needed); }

    std::vector<std::vector<int>> send_found(mpi::size());
    std::vector<std::vector<int>> recv_found(mpi::size());

    for (idx_t jpart = 0; jpart < nparts; ++jpart) {
        const std::vector<uid_t>& recv_node = recv_needed[proc[jpart]];
        const idx_t nb_recv_nodes           = idx_t(recv_node.size()) / varsize;
        for (idx_t jnode = 0; jnode < nb_recv_nodes; ++jnode) {
            uid_t uid = recv_node[jnode * varsize + 0];
            int inode = recv_node[jnode * varsize + 1];
            send_found[proc[jpart]].push_back(inode);
            send_found[proc[jpart]].push_back(lookup.count(uid) ? lookup[uid] : -1);
        }
    }

    ATLAS_TRACE_MPI(ALLTOALL) { mpi::comm().allToAll(send_found, recv_found); }

    std::stringstream errstream;
    size_t failed{0};
    const auto gidx = array::make_view<gidx_t, 1>(nodes.global_index());
    for (idx_t jpart = 0; jpart < nparts; ++jpart) {
        const std::vector<int>& recv_node = recv_found[proc[jpart]];
        const idx_t nb_recv_nodes         = recv_node.size() / 2;
        // array::ArrayView<int,2> recv_node( recv_found[ proc[jpart] ].data(),
        //     array::make_shape(recv_found[ proc[jpart] ].size()/2,2) );
        for (idx_t jnode = 0; jnode < nb_recv_nodes; ++jnode) {
            idx_t inode      = recv_node[jnode * 2 + 0];
            idx_t ridx_inode = recv_node[jnode * 2 + 1];
            if (ridx_inode >= 0) {
                ridx(recv_node[jnode * 2 + 0]) = ridx_inode;
            }
            else {
                ++failed;
                errstream << "\n[" << mpi::rank() << "] "
                          << "Node with global index " << gidx(inode) << " not found on part [" << part(inode) << "]";
            }
        }
    }

    if (failed) {
        throw_AssertionFailed(errstream.str(), Here());
    }

    return nodes.remote_index();
}

//----------------------------------------------------------------------------------------------------------------------

Field& build_nodes_partition(mesh::Nodes& nodes) {
    ATLAS_TRACE();
    return nodes.partition();
}

//----------------------------------------------------------------------------------------------------------------------

Field& build_edges_partition(Mesh& mesh) {
    ATLAS_TRACE();

    const mesh::Nodes& nodes = mesh.nodes();

    idx_t mypart = mpi::rank();

    mesh::HybridElements& edges = mesh.edges();
    auto edge_part              = array::make_view<int, 1>(edges.partition());
    //const auto edge_glb_idx     = array::make_view<gidx_t, 1>( edges.global_index() );

    const auto& edge_nodes   = edges.node_connectivity();
    const auto& edge_to_elem = edges.cell_connectivity();

    const auto node_part = array::make_view<int, 1>(nodes.partition());
    const auto xy        = array::make_view<double, 2>(nodes.xy());
    const auto flags     = array::make_view<int, 1>(nodes.flags());
    const auto node_gidx = array::make_view<gidx_t, 1>(nodes.global_index());

    const auto elem_gidx = array::make_view<gidx_t, 1>(mesh.cells().global_index());
    const auto elem_part = array::make_view<int, 1>(mesh.cells().partition());
    //const auto elem_halo = array::make_view<int, 1>( mesh.cells().halo() );

    auto check_flags = [&](idx_t jedge, int flag) {
        idx_t ip1 = edge_nodes(jedge, 0);
        idx_t ip2 = edge_nodes(jedge, 1);
        return Topology::check(flags(ip1), flag) && Topology::check(flags(ip2), flag);
    };
    auto domain_bdry = [&](idx_t jedge) {
        if (check_flags(jedge, Topology::BC | Topology::NORTH)) {
            return true;
        }
        if (check_flags(jedge, Topology::BC | Topology::SOUTH)) {
            return true;
        }
        if (check_flags(jedge, Topology::BC | Topology::WEST)) {
            return true;
        }
        if (check_flags(jedge, Topology::BC | Topology::EAST)) {
            return true;
        }
        return false;
    };
    auto periodic_east = [&](idx_t jedge) {
        if (check_flags(jedge, Topology::PERIODIC | Topology::EAST)) {
            return true;
        }
        return false;
    };
    auto periodic_west = [&](idx_t jedge) {
        if (check_flags(jedge, Topology::PERIODIC | Topology::WEST)) {
            return true;
        }
        return false;
    };

    auto periodic_west_bdry = [&](idx_t jedge) {
        if (check_flags(jedge, Topology::PERIODIC | Topology::WEST | Topology::BC)) {
            return true;
        }
        return false;
    };


    idx_t nb_edges = edges.size();

    std::vector<gidx_t> bdry_edges;
    bdry_edges.reserve(nb_edges);
    std::map<gidx_t, idx_t> global_to_local;


    PeriodicTransform transform_periodic_east(-360.);
    PeriodicTransform transform_periodic_west(+360.);
    UniqueLonLat _compute_uid(mesh);
    auto compute_uid = [&](idx_t jedge) -> gidx_t {
        if (periodic_east(jedge)) {
            return -_compute_uid(edge_nodes.row(jedge), transform_periodic_east);
        }
        else if (periodic_west_bdry(jedge)) {
            return _compute_uid(edge_nodes.row(jedge));
        }
        else if (periodic_west(jedge)) {
            return -_compute_uid(edge_nodes.row(jedge), transform_periodic_west);
        }
        else {
            return _compute_uid(edge_nodes.row(jedge));
        }
    };


    // should be unit-test
    {
        ATLAS_ASSERT(util::unique_lonlat(360., 0., transform_periodic_east) == util::unique_lonlat(0., 0.));
        ATLAS_ASSERT(util::unique_lonlat(0., 0., transform_periodic_west) == util::unique_lonlat(360., 0.));
    }


    for (idx_t jedge = 0; jedge < nb_edges; ++jedge) {
        gidx_t edge_gidx           = compute_uid(jedge);
        global_to_local[edge_gidx] = jedge;

#ifdef DEBUGGING_PARFIELDS
        if (FIND_EDGE(jedge)) {
            std::cout << "[" << mypart << "] " << EDGE(jedge) << " has gidx " << edge_gidx << std::endl;
            if (periodic_east(jedge)) {
                std::cout << "[" << mypart << "] " << EDGE(jedge) << " is periodic_east " << std::endl;
            }
            else if (periodic_west(jedge)) {
                std::cout << "[" << mypart << "] " << EDGE(jedge) << " is periodic_west " << std::endl;
            }
            else {
                std::cout << "[" << mypart << "] " << EDGE(jedge) << " is not periodic" << std::endl;
            }
        }
        if (FIND_GIDX(edge_gidx))
            std::cout << "[" << mypart << "] "
                      << "has " << EDGE(jedge) << " with gidx " << edge_gidx << std::endl;
#endif

        idx_t ip1 = edge_nodes(jedge, 0);
        idx_t ip2 = edge_nodes(jedge, 1);
        int pn1   = node_part(ip1);
        int pn2   = node_part(ip2);
        int y1    = util::microdeg(xy(ip1, YY));
        int y2    = util::microdeg(xy(ip2, YY));
        int p;
        if (y1 == y2) {
            int x1 = util::microdeg(xy(ip1, XX));
            int x2 = util::microdeg(xy(ip2, XX));
            p      = (x1 < x2) ? pn1 : pn2;
        }
        else {
            p = (y1 > y2) ? pn1 : pn2;
        }

        idx_t elem1   = edge_to_elem(jedge, 0);
        idx_t elem2   = edge_to_elem(jedge, 1);
        idx_t missing = edge_to_elem.missing_value();
        if (elem1 == missing && elem2 == missing) {
            // Don't attempt to change p
        }
        else if (elem1 == missing) {
            ATLAS_NOTIMPLEMENTED;
        }
        else if (elem2 == missing) {
            if (pn1 == pn2) {
                p = pn1;
            }
            else if (periodic_east(jedge)) {
#ifdef DEBUGGING_PARFIELDS
                if (FIND_EDGE(jedge))
                    std::cout << "[" << mypart << "] "
                              << "periodic_east" << std::endl;
#endif
                bdry_edges.push_back(edge_gidx);
                p = -1;
            }
            else if (periodic_west_bdry(jedge)) {
#ifdef DEBUGGING_PARFIELDS
                if (FIND_EDGE(jedge))
                    std::cout << "[" << mypart << "] "
                              << "periodic_west_bdry" << std::endl;
#endif
                p = elem_part(elem1);
            }
            else if (periodic_west(jedge)) {
#ifdef DEBUGGING_PARFIELDS
                if (FIND_EDGE(jedge))
                    std::cout << "[" << mypart << "] "
                              << "periodic_west" << std::endl;
#endif
                bdry_edges.push_back(edge_gidx);
                p = -1;
            }
            else if (domain_bdry(jedge)) {
#ifdef DEBUGGING_PARFIELDS
                if (FIND_EDGE(jedge))
                    std::cout << "[" << mypart << "] "
                              << "domain_bdry" << std::endl;
#endif
                p = elem_part(elem1);
            }
            else {
                bdry_edges.push_back(edge_gidx);
                p = -1;
            }
        }
        else if (p != elem_part(elem1) && p != elem_part(elem2)) {
            p = (p == pn1) ? pn2 : pn1;

            if (p != elem_part(elem1) and p != elem_part(elem2)) {
                std::stringstream msg;
                msg << "[" << mpi::rank() << "] " << EDGE(jedge) << " has nodes and elements of different rank: elem1[p"
                    << elem_part(elem1) << "] elem2[p" << elem_part(elem2) << "]";
                throw_Exception(msg.str(), Here());
            }
        }
        edge_part(jedge) = p;
    }
    std::sort(bdry_edges.begin(), bdry_edges.end());
    auto is_bdry_edge = [&bdry_edges](gidx_t gid) {
        std::vector<uid_t>::iterator it = std::lower_bound(bdry_edges.begin(), bdry_edges.end(), gid);
        bool found                      = !(it == bdry_edges.end() || gid < *it);
        return found;
    };

    int mpi_size = mpi::size();
    mpi::Buffer<gidx_t, 1> recv_bdry_edges_from_parts(mpi_size);
    std::vector<std::vector<gidx_t>> send_gidx(mpi_size);
    std::vector<std::vector<int>> send_part(mpi_size);
    std::vector<std::vector<gidx_t>> send_bdry_gidx(mpi_size);
    std::vector<std::vector<int>> send_bdry_elem_part(mpi_size);
    std::vector<std::vector<gidx_t>> send_bdry_elem_gidx(mpi_size);
    mpi::comm().allGatherv(bdry_edges.begin(), bdry_edges.end(), recv_bdry_edges_from_parts);
    for (int p = 0; p < mpi_size; ++p) {
        auto view = recv_bdry_edges_from_parts[p];
        for (size_t j = 0; j < view.size(); ++j) {
            gidx_t gidx        = view[j];
            gidx_t master_gidx = std::abs(gidx);
            if (global_to_local.count(master_gidx)) {
                idx_t iedge = global_to_local[master_gidx];
#ifdef DEBUGGING_PARFIELDS
                if (FIND_GIDX(master_gidx))
                    std::cout << "[" << mypart << "] found " << EDGE(iedge) << std::endl;
#endif
                if (not is_bdry_edge(master_gidx)) {
                    send_gidx[p].push_back(gidx);
                    send_part[p].push_back(edge_part(iedge));
#ifdef DEBUGGING_PARFIELDS
                    if (FIND_EDGE(iedge)) {
                        std::cout << "[" << mypart << "] found " << EDGE(iedge) " for part " << p << std::endl;
                    }
#endif
                }
                else {  // boundary edge with nodes of different rank
                    idx_t ielem = (edge_to_elem(iedge, 0) != edge_to_elem.missing_value() ? edge_to_elem(iedge, 0)
                                                                                          : edge_to_elem(iedge, 1));
                    send_bdry_gidx[p].push_back(gidx);
                    send_bdry_elem_part[p].push_back(elem_part(ielem));
                    send_bdry_elem_gidx[p].push_back(elem_gidx(ielem));
                }
            }
        }
    }
    std::vector<std::vector<gidx_t>> recv_gidx(mpi_size);
    std::vector<std::vector<int>> recv_part(mpi_size);

    mpi::comm().allToAll(send_gidx, recv_gidx);
    mpi::comm().allToAll(send_part, recv_part);
    for (int p = 0; p < mpi_size; ++p) {
        const auto& recv_gidx_p = recv_gidx[p];
        const auto& recv_part_p = recv_part[p];
        for (size_t j = 0; j < recv_gidx_p.size(); ++j) {
            idx_t iedge = global_to_local[recv_gidx_p[j]];
            // int prev           = edge_part( iedge );
            edge_part(iedge) = recv_part_p[j];
            // if( edge_part(iedge) != prev )
            //   Log::error() << EDGE(iedge) << " part " << prev << " --> " <<
            //   edge_part(iedge) << std::endl;
        }
    }

    std::vector<std::vector<gidx_t>> recv_bdry_gidx(mpi_size);
    std::vector<std::vector<int>> recv_bdry_elem_part(mpi_size);
    std::vector<std::vector<gidx_t>> recv_bdry_elem_gidx(mpi_size);
    mpi::comm().allToAll(send_bdry_gidx, recv_bdry_gidx);
    mpi::comm().allToAll(send_bdry_elem_part, recv_bdry_elem_part);
    mpi::comm().allToAll(send_bdry_elem_gidx, recv_bdry_elem_gidx);
    for (int p = 0; p < mpi_size; ++p) {
        const auto& recv_bdry_gidx_p      = recv_bdry_gidx[p];
        const auto& recv_bdry_elem_part_p = recv_bdry_elem_part[p];
        const auto& recv_bdry_elem_gidx_p = recv_bdry_elem_gidx[p];
        for (size_t j = 0; j < recv_bdry_gidx_p.size(); ++j) {
            idx_t iedge = global_to_local[recv_bdry_gidx_p[j]];
            idx_t e1    = (edge_to_elem(iedge, 0) != edge_to_elem.missing_value() ? edge_to_elem(iedge, 0)
                                                                                  : edge_to_elem(iedge, 1));
            if (elem_gidx(e1) != recv_bdry_elem_gidx_p[j]) {
                idx_t ip1 = edge_nodes(iedge, 0);
                idx_t ip2 = edge_nodes(iedge, 1);
                int pn1   = node_part(ip1);
                int pn2   = node_part(ip2);
                int y1    = util::microdeg(xy(ip1, YY));
                int y2    = util::microdeg(xy(ip2, YY));
                int ped;
                if (y1 == y2) {
                    int x1 = util::microdeg(xy(ip1, XX));
                    int x2 = util::microdeg(xy(ip2, XX));
                    ped    = (x1 < x2) ? pn1 : pn2;
                }
                else {
                    ped = (y1 > y2) ? pn1 : pn2;
                }
                int pe1 = elem_part(e1);
                int pe2 = recv_bdry_elem_part_p[j];
                if (ped != pe1 && ped != pe2) {
                    ped = (ped == pn1) ? pn2 : pn1;
                    if (ped != pe1 && p != pe2) {
                        std::stringstream msg;
                        msg << "[" << mpi::rank() << "] " << EDGE(iedge)
                            << " has nodes and elements of different rank: elem1[p" << pe1 << "] elem2[p" << pe2 << "]";
                        throw_Exception(msg.str(), Here());
                    }
                }
                edge_part(iedge) = ped;
            }
        }
    }

    // Sanity check
    auto edge_flags     = array::make_view<int, 1>(edges.flags());
    auto is_pole_edge   = [&](idx_t e) { return Topology::check(edge_flags(e), Topology::POLE); };
    bool has_pole_edges = false;
    mesh.edges().metadata().get("pole_edges", has_pole_edges);
    int insane = 0;
    for (idx_t jedge = 0; jedge < nb_edges; ++jedge) {
        idx_t ip1   = edge_nodes(jedge, 0);
        idx_t ip2   = edge_nodes(jedge, 1);
        idx_t elem1 = edge_to_elem(jedge, 0);
        idx_t elem2 = edge_to_elem(jedge, 1);
        int p       = edge_part(jedge);
        int pn1     = node_part(ip1);
        int pn2     = node_part(ip2);
        if (has_pole_edges && is_pole_edge(jedge)) {
            if (p != pn1 || p != pn2) {
                Log::error() << "pole edge " << EDGE(jedge) << " [p" << p << "] is not correct" << std::endl;
                insane = 1;
            }
        }
        else {
            if (elem1 == edge_to_elem.missing_value() && elem2 == edge_to_elem.missing_value()) {
                Log::error() << EDGE(jedge) << " has no neighbouring elements" << std::endl;
                insane = 1;
            }
        }
        bool edge_is_partition_boundary =
            (elem1 == edge_to_elem.missing_value() || elem2 == edge_to_elem.missing_value());
        bool edge_partition_is_same_as_one_of_nodes = (p == pn1 || p == pn2);
        if (edge_is_partition_boundary) {
            if (not edge_partition_is_same_as_one_of_nodes) {
                gidx_t edge_gidx = compute_uid(jedge);

                if (elem1 != edge_to_elem.missing_value()) {
                    Log::error() << "[" << mypart << "] " << EDGE(jedge) << " " << edge_gidx << " [p" << p
                                 << "] at partition_boundary is not correct. elem1[p" << elem_part(elem1) << "]"
                                 << std::endl;
                }
                else {
                    Log::error() << "[" << mypart << "] " << EDGE(jedge) << " " << edge_gidx << " [p" << p
                                 << "] at partition_boundary is not correct elem2[p" << elem_part(elem2) << "]"
                                 << std::endl;
                }
                insane = 1;
            }
        }
        else {
            int pe1                                     = elem_part(elem1);
            int pe2                                     = elem_part(elem2);
            bool edge_partition_is_same_as_one_of_elems = (p == pe1 || p == pe2);
            if (not edge_partition_is_same_as_one_of_elems and not edge_partition_is_same_as_one_of_nodes) {
                Log::error() << EDGE(jedge) << " is not correct elem1[p" << pe1 << "] elem2[p" << pe2 << "]"
                             << std::endl;
                insane = 1;
            }
        }
    }
    mpi::comm().allReduceInPlace(insane, eckit::mpi::max());
    if (insane && mpi::rank() == 0) {
        throw_Exception("Sanity check failed", Here());
    }

    //#ifdef DEBUGGING_PARFIELDS
    //        if( OWNED_EDGE(jedge) )
    //          DEBUG( EDGE(jedge) << " -->    part " << edge_part(jedge));
    //#endif

    //#ifdef DEBUGGING_PARFIELDS_disable
    //    if( PERIODIC_EDGE(jedge) )
    //      DEBUG_VAR( "           the part is " << edge_part(jedge) );
    //#endif
    //  }

    return edges.partition();
}

Field& build_edges_remote_idx(Mesh& mesh) {
    ATLAS_TRACE();

    const mesh::Nodes& nodes = mesh.nodes();
    UniqueLonLat compute_uid(mesh);

    idx_t mypart = mpi::rank();
    idx_t nparts = mpi::size();

    mesh::HybridElements& edges = mesh.edges();

    auto edge_ridx = array::make_indexview<idx_t, 1>(edges.remote_index());

    const array::ArrayView<int, 1> edge_part             = array::make_view<int, 1>(edges.partition());
    const mesh::HybridElements::Connectivity& edge_nodes = edges.node_connectivity();

    array::ArrayView<const double, 2> xy = array::make_view<double, 2>(nodes.xy());
    array::ArrayView<const int, 1> flags = array::make_view<int, 1>(nodes.flags());
#ifdef DEBUGGING_PARFIELDS
    array::ArrayView<gidx_t, 1> node_gidx = array::make_view<gidx_t, 1>(nodes.global_index());
    array::ArrayView<int, 1> node_part    = array::make_view<int, 1>(nodes.partition());
#endif

    auto edge_flags     = array::make_view<int, 1>(edges.flags());
    auto is_pole_edge   = [&](idx_t e) { return Topology::check(edge_flags(e), Topology::POLE); };
    bool has_pole_edges = false;
    mesh.edges().metadata().get("pole_edges", has_pole_edges);

    const int nb_edges = edges.size();

    double centroid[2];
    std::vector<std::vector<uid_t>> send_needed(mpi::size());
    std::vector<std::vector<uid_t>> recv_needed(mpi::size());
    int sendcnt = 0;
    std::map<uid_t, int> lookup;

    PeriodicTransform transform;

    for (int jedge = 0; jedge < nb_edges; ++jedge) {
        int ip1      = edge_nodes(jedge, 0);
        int ip2      = edge_nodes(jedge, 1);
        centroid[XX] = 0.5 * (xy(ip1, XX) + xy(ip2, XX));
        centroid[YY] = 0.5 * (xy(ip1, YY) + xy(ip2, YY));
        if (has_pole_edges && is_pole_edge(jedge)) {
            centroid[YY] = centroid[YY] > 0 ? 90. : -90.;
        }

        bool needed(false);

        if ((Topology::check(flags(ip1), Topology::PERIODIC) &&
             !Topology::check(flags(ip1), Topology::BC | Topology::WEST) &&
             Topology::check(flags(ip2), Topology::PERIODIC) &&
             !Topology::check(flags(ip2), Topology::BC | Topology::WEST)) ||
            (Topology::check(flags(ip1), Topology::PERIODIC) &&
             !Topology::check(flags(ip1), Topology::BC | Topology::WEST) &&
             Topology::check(flags(ip2), Topology::BC | Topology::WEST)) ||
            (Topology::check(flags(ip1), Topology::BC | Topology::WEST) &&
             Topology::check(flags(ip2), Topology::PERIODIC) &&
             !Topology::check(flags(ip2), Topology::BC | Topology::WEST))) {
            needed = true;
            if (Topology::check(flags(ip1), Topology::EAST)) {
                transform(centroid, -1);
            }
            else {
                transform(centroid, +1);
            }
        }

        uid_t uid = util::unique_lonlat(centroid);
        if (idx_t(edge_part(jedge)) == mypart && !needed)  // All interior edges fall here
        {
            lookup[uid]      = jedge;
            edge_ridx(jedge) = jedge;

#ifdef DEBUGGING_PARFIELDS
            if (FIND_EDGE(jedge)) {
                ATLAS_DEBUG("Found " << EDGE(jedge));
            }
#endif
        }
        else  // All ghost edges PLUS the periodic edges identified edges above
              // fall here
        {
            send_needed[edge_part(jedge)].push_back(uid);
            send_needed[edge_part(jedge)].push_back(jedge);
#ifdef DEBUGGING_PARFIELDS
            send_needed[edge_part(jedge)].push_back(node_gidx(ip1));
            send_needed[edge_part(jedge)].push_back(node_gidx(ip2));
            send_needed[edge_part(jedge)].push_back(node_part(ip1));
            send_needed[edge_part(jedge)].push_back(node_part(ip2));
#endif
            sendcnt++;
        }
    }

    idx_t varsize = 2;
#ifdef DEBUGGING_PARFIELDS
    varsize = 6;
#endif

    ATLAS_TRACE_MPI(ALLTOALL) { mpi::comm().allToAll(send_needed, recv_needed); }

    std::vector<std::vector<int>> send_found(mpi::size());
    std::vector<std::vector<int>> recv_found(mpi::size());

    std::map<uid_t, int>::iterator found;
    for (idx_t jpart = 0; jpart < nparts; ++jpart) {
        const std::vector<uid_t>& recv_edge = recv_needed[jpart];
        const idx_t nb_recv_edges           = idx_t(recv_edge.size()) / varsize;
        // array::ArrayView<uid_t,2> recv_edge( recv_needed[ jpart ].data(),
        //     array::make_shape(recv_needed[ jpart ].size()/varsize,varsize) );
        for (idx_t jedge = 0; jedge < nb_recv_edges; ++jedge) {
            uid_t recv_uid = recv_edge[jedge * varsize + 0];
            int recv_idx   = recv_edge[jedge * varsize + 1];
            found          = lookup.find(recv_uid);
            if (found != lookup.end()) {
                send_found[jpart].push_back(recv_idx);
                send_found[jpart].push_back(found->second);
            }
            else {
                std::stringstream msg;
#ifdef DEBUGGING_PARFIELDS
                msg << "Edge(" << recv_edge[jedge * varsize + 2] << "[p" << recv_edge[jedge * varsize + 4] << "] "
                    << recv_edge[jedge * varsize + 3] << "[p" << recv_edge[jedge * varsize + 5] << "])";
#else
                msg << "Edge with uid " << recv_uid;
#endif
                msg << " requested by rank [" << jpart << "]";
                msg << " that should be owned by " << mpi::rank()
                    << " is not found. This could be because no "
                       "halo was built.";
                // throw_Exception(msg.str(),Here());
                Log::warning() << msg.str() << " @ " << Here() << std::endl;
            }
        }
    }

    ATLAS_TRACE_MPI(ALLTOALL) { mpi::comm().allToAll(send_found, recv_found); }

    for (idx_t jpart = 0; jpart < nparts; ++jpart) {
        const std::vector<int>& recv_edge = recv_found[jpart];
        const idx_t nb_recv_edges         = recv_edge.size() / 2;
        // array::ArrayView<int,2> recv_edge( recv_found[ jpart ].data(),
        //     array::make_shape(recv_found[ jpart ].size()/2,2) );
        for (idx_t jedge = 0; jedge < nb_recv_edges; ++jedge) {
            edge_ridx(recv_edge[jedge * 2 + 0]) = recv_edge[jedge * 2 + 1];
        }
    }
    return edges.remote_index();
}

Field& build_edges_global_idx(Mesh& mesh) {
    ATLAS_TRACE();

    UniqueLonLat compute_uid(mesh);

    int nparts = mpi::size();
    idx_t root = 0;

    mesh::HybridElements& edges = mesh.edges();

    array::make_view<gidx_t, 1>(edges.global_index()).assign(-1);
    array::ArrayView<gidx_t, 1> edge_gidx = array::make_view<gidx_t, 1>(edges.global_index());

    const mesh::HybridElements::Connectivity& edge_nodes = edges.node_connectivity();
    array::ArrayView<double, 2> xy                       = array::make_view<double, 2>(mesh.nodes().xy());
    auto edge_flags                                      = array::make_view<int, 1>(edges.flags());
    auto is_pole_edge   = [&](idx_t e) { return Topology::check(edge_flags(e), Topology::POLE); };
    bool has_pole_edges = false;
    mesh.edges().metadata().get("pole_edges", has_pole_edges);

    /*
 * Sorting following edge_gidx will define global order of
 * gathered fields. Special care needs to be taken for
 * pole edges, as their centroid might coincide with
 * other edges
 */
    double centroid[2];
    int nb_edges = edges.size();
    for (int jedge = 0; jedge < nb_edges; ++jedge) {
        if (edge_gidx(jedge) <= 0) {
            centroid[XX] = 0.5 * (xy(edge_nodes(jedge, 0), XX) + xy(edge_nodes(jedge, 1), XX));
            centroid[YY] = 0.5 * (xy(edge_nodes(jedge, 0), YY) + xy(edge_nodes(jedge, 1), YY));
            if (has_pole_edges && is_pole_edge(jedge)) {
                centroid[YY] = centroid[YY] > 0 ? 90. : -90.;
            }
            edge_gidx(jedge) = util::unique_lonlat(centroid);
        }
    }

    /*
 * REMOTE INDEX BASE = 1
 */

    // 1) Gather all global indices, together with location
    array::ArrayT<uid_t> loc_edge_id_arr(nb_edges);
    array::ArrayView<uid_t, 1> loc_edge_id = array::make_view<uid_t, 1>(loc_edge_id_arr);

    for (int jedge = 0; jedge < nb_edges; ++jedge) {
        loc_edge_id(jedge) = edge_gidx(jedge);
    }

    std::vector<int> recvcounts(mpi::size());
    std::vector<int> recvdispls(mpi::size());

    ATLAS_TRACE_MPI(GATHER) { mpi::comm().gather(nb_edges, recvcounts, root); }

    recvdispls[0] = 0;
    for (int jpart = 1; jpart < nparts; ++jpart)  // start at 1
    {
        recvdispls[jpart] = recvcounts[jpart - 1] + recvdispls[jpart - 1];
    }
    int glb_nb_edges = std::accumulate(recvcounts.begin(), recvcounts.end(), 0);

    array::ArrayT<uid_t> glb_edge_id_arr(glb_nb_edges);
    array::ArrayView<uid_t, 1> glb_edge_id = array::make_view<uid_t, 1>(glb_edge_id_arr);

    ATLAS_TRACE_MPI(GATHER) {
        mpi::comm().gatherv(loc_edge_id.data(), loc_edge_id.size(), glb_edge_id.data(), recvcounts.data(),
                            recvdispls.data(), root);
    }

    // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
    std::vector<Node> edge_sort;
    edge_sort.reserve(glb_nb_edges);
    for (idx_t jedge = 0; jedge < glb_edge_id.shape(0); ++jedge) {
        edge_sort.emplace_back(glb_edge_id(jedge), jedge);
    }
    std::sort(edge_sort.begin(), edge_sort.end());

    // Assume edge gid start
    uid_t gid(0);
    for (size_t jedge = 0; jedge < edge_sort.size(); ++jedge) {
        if (jedge == 0) {
            ++gid;
        }
        else if (edge_sort[jedge].g != edge_sort[jedge - 1].g) {
            ++gid;
        }
        int iedge          = edge_sort[jedge].i;
        glb_edge_id(iedge) = gid;
    }

    // 3) Scatter renumbered back
    ATLAS_TRACE_MPI(SCATTER) {
        mpi::comm().scatterv(glb_edge_id.data(), recvcounts.data(), recvdispls.data(), loc_edge_id.data(),
                             loc_edge_id.size(), root);
    }

    for (int jedge = 0; jedge < nb_edges; ++jedge) {
        edge_gidx(jedge) = loc_edge_id(jedge);
    }

    return edges.global_index();
}

//----------------------------------------------------------------------------------------------------------------------

Field& build_cells_remote_idx(mesh::Cells& cells, const mesh::Nodes& nodes) {
    ATLAS_TRACE();
    idx_t mypart = static_cast<idx_t>(mpi::rank());
    idx_t nparts = static_cast<idx_t>(mpi::size());

    // This piece should be somewhere central ... could be NPROMA ?
    // ---------->
    std::vector<idx_t> proc(nparts);
    for (idx_t jpart = 0; jpart < nparts; ++jpart) {
        proc[jpart] = jpart;
    }
    // <---------

    auto ridx                 = array::make_indexview<idx_t, 1>(cells.remote_index());
    const auto part           = array::make_view<int, 1>(cells.partition());
    const auto gidx           = array::make_view<gidx_t, 1>(cells.global_index());
    const auto flags          = array::make_view<int, 1>(cells.flags());
    const auto& element_nodes = cells.node_connectivity();
    idx_t nb_cells            = cells.size();

    const PeriodicTransform transform_periodic_east(-360.);
    const PeriodicTransform transform_periodic_west(+360.);
    const UniqueLonLat compute_uid_lonlat(nodes);
    auto compute_uid = [&](idx_t jcell) {
        constexpr int PERIODIC = util::Topology::PERIODIC;
        constexpr int EAST     = util::Topology::EAST;
        constexpr int WEST     = util::Topology::WEST;
        const auto flags_view  = util::Bitflags::view(flags(jcell));
        if (flags_view.check(PERIODIC | EAST)) {
            return compute_uid_lonlat(element_nodes.row(jcell), transform_periodic_east);
        }
        if (flags_view.check(PERIODIC | WEST)) {
            return compute_uid_lonlat(element_nodes.row(jcell), transform_periodic_west);
        }
        return compute_uid_lonlat(element_nodes.row(jcell));
    };

    idx_t varsize = 2;

    std::vector<std::vector<uid_t>> send_needed(mpi::size());
    std::vector<std::vector<uid_t>> recv_needed(mpi::size());
    int sendcnt = 0;
    std::map<uid_t, int> lookup;
    for (idx_t jcell = 0; jcell < nb_cells; ++jcell) {
        uid_t uid = compute_uid(jcell);

        if (idx_t(part(jcell)) == mypart) {
            lookup[uid] = jcell;
            ridx(jcell) = jcell;
        }
        else {
            ATLAS_ASSERT(jcell < part.shape(0));
            if (part(jcell) >= static_cast<int>(proc.size())) {
                std::stringstream msg;
                msg << "Assertion [part(" << jcell << ") < proc.size()] failed\n"
                    << "part(" << jcell << ") = " << part(jcell) << "\n"
                    << "proc.size() = " << proc.size();
                throw_AssertionFailed(msg.str(), Here());
            }
            ATLAS_ASSERT(part(jcell) < (idx_t)proc.size());
            ATLAS_ASSERT((size_t)proc[part(jcell)] < send_needed.size());
            send_needed[proc[part(jcell)]].push_back(uid);
            send_needed[proc[part(jcell)]].push_back(jcell);
            sendcnt++;
        }
    }

    ATLAS_TRACE_MPI(ALLTOALL) { mpi::comm().allToAll(send_needed, recv_needed); }

    std::vector<std::vector<int>> send_found(mpi::size());
    std::vector<std::vector<int>> recv_found(mpi::size());

    for (idx_t jpart = 0; jpart < nparts; ++jpart) {
        const std::vector<uid_t>& recv_cell = recv_needed[proc[jpart]];
        const idx_t nb_recv_cells           = idx_t(recv_cell.size()) / varsize;
        for (idx_t jcell = 0; jcell < nb_recv_cells; ++jcell) {
            uid_t uid = recv_cell[jcell * varsize + 0];
            int icell = recv_cell[jcell * varsize + 1];
            send_found[proc[jpart]].push_back(icell);
            send_found[proc[jpart]].push_back(lookup.count(uid) ? lookup[uid] : -1);
        }
    }

    ATLAS_TRACE_MPI(ALLTOALL) { mpi::comm().allToAll(send_found, recv_found); }

    std::stringstream errstream;
    size_t failed{0};
    for (idx_t jpart = 0; jpart < nparts; ++jpart) {
        const std::vector<int>& recv_cell = recv_found[proc[jpart]];
        const idx_t nb_recv_cells         = recv_cell.size() / 2;
        // array::ArrayView<int,2> recv_node( recv_found[ proc[jpart] ].data(),
        //     array::make_shape(recv_found[ proc[jpart] ].size()/2,2) );
        for (idx_t jcell = 0; jcell < nb_recv_cells; ++jcell) {
            idx_t icell      = recv_cell[jcell * 2 + 0];
            idx_t ridx_icell = recv_cell[jcell * 2 + 1];
            if (ridx_icell >= 0) {
                ridx(icell) = ridx_icell;
            }
            else {
                ++failed;
                errstream << "\n[" << mpi::rank() << "] "
                          << "Cell " << gidx(icell) << " not found on part [" << part(icell) << "]";
            }
        }
    }
    if (failed) {
        throw_AssertionFailed(errstream.str(), Here());
    }
    return cells.remote_index();
}

void build_cells_parallel_fields(Mesh& mesh) {
    mpi::Scope mpi_scope(mesh.mpi_comm());

    bool parallel = false;
    mesh.cells().metadata().get("parallel", parallel);
    if (!parallel) {
        build_cells_remote_idx(mesh.cells(), mesh.nodes());
    }

    mesh.cells().metadata().set("parallel", true);
}

//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_parallel_fields(Mesh::Implementation* mesh) {
    ATLAS_ASSERT(mesh != nullptr, "Cannot access uninitialised atlas_Mesh");
    Mesh m(mesh);
    build_parallel_fields(m);
}
void atlas__build_nodes_parallel_fields(mesh::Nodes* nodes) {
    ATLAS_ASSERT(nodes != nullptr, "Cannot access uninitialised atlas_mesh_Nodes");
    build_nodes_parallel_fields(*nodes);
}

void atlas__build_edges_parallel_fields(Mesh::Implementation* mesh) {
    ATLAS_ASSERT(mesh != nullptr, "Cannot access uninitialised atlas_Mesh");
    Mesh m(mesh);
    build_edges_parallel_fields(m);
}

void atlas__renumber_nodes_glb_idx(mesh::Nodes* nodes) {
    ATLAS_ASSERT(nodes != nullptr, "Cannot access uninitialised atlas_mesh_Nodes");
    renumber_nodes_glb_idx(*nodes);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
