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
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "atlas/array.h"
#include "atlas/array/IndexView.h"
#include "atlas/field/Field.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/detail/AccumulateFacets.h"
#include "atlas/parallel/mpi/Buffer.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/parallel/omp/sort.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/MicroDeg.h"
#include "atlas/util/PeriodicTransform.h"
#include "atlas/util/Unique.h"

//#define DEBUG_OUTPUT
#ifdef DEBUG_OUTPUT
#include "atlas/mesh/actions/BuildParallelFields.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/output/Gmsh.h"
#endif

// #define ATLAS_103
// #define ATLAS_103_SORT

using atlas::mesh::detail::accumulate_facets;
using atlas::util::LonLatMicroDeg;
using atlas::util::microdeg;
using atlas::util::PeriodicTransform;
using atlas::util::UniqueLonLat;
using Topology = atlas::mesh::Nodes::Topology;

namespace atlas {
namespace mesh {
namespace actions {

struct Entity {
    Entity(gidx_t gid, idx_t idx) {
        g = gid;
        i = idx;
    }
    gidx_t g;
    idx_t i;
    bool operator<(const Entity& other) const { return (g < other.g); }
};

void make_nodes_global_index_human_readable(const mesh::actions::BuildHalo& build_halo, mesh::Nodes& nodes,
                                            bool do_all) {
    ATLAS_TRACE();
    // TODO: ATLAS-14: fix renumbering of EAST periodic boundary points
    // --> Those specific periodic points at the EAST boundary are not checked for
    // uid,
    //     and could receive different gidx for different tasks

    // unused // int mypart = mpi::rank();
    int nparts  = mpi::size();
    size_t root = 0;

    array::ArrayView<gidx_t, 1> nodes_glb_idx = array::make_view<gidx_t, 1>(nodes.global_index());
    // nodes_glb_idx.dump( Log::info() );
    //  ATLAS_DEBUG( "min = " << nodes.global_index().metadata().getLong("min") );
    //  ATLAS_DEBUG( "max = " << nodes.global_index().metadata().getLong("max") );
    //  ATLAS_DEBUG( "human_readable = " <<
    //  nodes.global_index().metadata().getBool("human_readable") );
    gidx_t glb_idx_max = 0;

    std::vector<int> points_to_edit;

    if (do_all) {
        points_to_edit.resize(nodes_glb_idx.size());
        for (idx_t i = 0; i < nodes_glb_idx.size(); ++i) {
            points_to_edit[i] = i;
        }
    }
    else {
        glb_idx_max = nodes.global_index().metadata().getLong("max", 0);
        points_to_edit.resize(build_halo.periodic_points_local_index_.size());
        for (size_t i = 0; i < points_to_edit.size(); ++i) {
            points_to_edit[i] = build_halo.periodic_points_local_index_[i];
        }
    }

    std::vector<gidx_t> glb_idx(points_to_edit.size());
    int nb_nodes = static_cast<int>(glb_idx.size());
    for (idx_t i = 0; i < nb_nodes; ++i) {
        glb_idx[i] = nodes_glb_idx(points_to_edit[i]);
    }

    //  ATLAS_DEBUG_VAR( points_to_edit );
    //  ATLAS_DEBUG_VAR( points_to_edit.size() );

    /*
 * Sorting following gidx will define global order of
 * gathered fields. Special care needs to be taken for
 * pole edges, as their centroid might coincide with
 * other edges
 */
    // Try to recover
    //  {
    //    UniqueLonLat compute_uid(nodes);
    //    for( int jnode=0; jnode<nb_nodes; ++jnode ) {
    //      if( glb_idx(jnode) <= 0 ) {
    //        glb_idx(jnode) = compute_uid(jnode);
    //      }
    //    }
    //  }

    // 1) Gather all global indices, together with location

    std::vector<int> recvcounts(mpi::size());
    std::vector<int> recvdispls(mpi::size());

    ATLAS_TRACE_MPI(GATHER) { mpi::comm().gather(nb_nodes, recvcounts, root); }
    idx_t glb_nb_nodes = std::accumulate(recvcounts.begin(), recvcounts.end(), 0);

    recvdispls[0] = 0;
    for (int jpart = 1; jpart < nparts; ++jpart) {  // start at 1
        recvdispls[jpart] = recvcounts[jpart - 1] + recvdispls[jpart - 1];
    }

    std::vector<gidx_t> glb_idx_gathered(glb_nb_nodes);
    ATLAS_TRACE_MPI(GATHER) {
        mpi::comm().gatherv(glb_idx.data(), glb_idx.size(), glb_idx_gathered.data(), recvcounts.data(),
                            recvdispls.data(), root);
    }

    // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
    std::vector<Entity> node_sort;
    node_sort.reserve(glb_nb_nodes);
    const idx_t nb_glb_idx_gathered = static_cast<idx_t>(glb_idx_gathered.size());
    for (idx_t jnode = 0; jnode < nb_glb_idx_gathered; ++jnode) {
        node_sort.emplace_back(glb_idx_gathered[jnode], jnode);
    }

    ATLAS_TRACE_SCOPE("sort on rank 0") { std::sort(node_sort.begin(), node_sort.end()); }

    gidx_t gid               = glb_idx_max + 1;
    const idx_t nb_node_sort = static_cast<idx_t>(node_sort.size());
    for (idx_t jnode = 0; jnode < nb_node_sort; ++jnode) {
        if (jnode > 0 && node_sort[jnode].g != node_sort[jnode - 1].g) {
            ++gid;
        }
        int inode               = node_sort[jnode].i;
        glb_idx_gathered[inode] = gid;
    }

    // 3) Scatter renumbered back
    ATLAS_TRACE_MPI(SCATTER) {
        mpi::comm().scatterv(glb_idx_gathered.data(), recvcounts.data(), recvdispls.data(), glb_idx.data(),
                             static_cast<int>(glb_idx.size()), root);
    }

    for (int jnode = 0; jnode < nb_nodes; ++jnode) {
        nodes_glb_idx(points_to_edit[jnode]) = glb_idx[jnode];
    }

    // nodes_glb_idx.dump( Log::info() );
    // Log::info() << std::endl;
    nodes.global_index().metadata().set("human_readable", true);
}

// -------------------------------------------------------------------------------------

void make_cells_global_index_human_readable(const mesh::actions::BuildHalo& build_halo, mesh::Cells& cells,
                                            bool do_all) {
    ATLAS_TRACE();

    int nparts = static_cast<idx_t>(mpi::size());
    idx_t root = 0;

    array::ArrayView<gidx_t, 1> cells_glb_idx = array::make_view<gidx_t, 1>(cells.global_index());
    //  ATLAS_DEBUG( "min = " << cells.global_index().metadata().getLong("min") );
    //  ATLAS_DEBUG( "max = " << cells.global_index().metadata().getLong("max") );
    //  ATLAS_DEBUG( "human_readable = " <<
    //  cells.global_index().metadata().getBool("human_readable") );
    gidx_t glb_idx_max = 0;

    std::vector<int> cells_to_edit;

    if (do_all) {
        cells_to_edit.resize(cells_glb_idx.size());
        for (idx_t i = 0; i < cells_glb_idx.size(); ++i) {
            cells_to_edit[i] = i;
        }
    }
    else {
        idx_t nb_cells_to_edit(0);
        for (const auto& new_cells : build_halo.periodic_cells_local_index_) {
            nb_cells_to_edit += new_cells.size();
        }
        cells_to_edit.resize(nb_cells_to_edit);
        int c{0};
        int i{0};
        for (int t = 0; t < cells.nb_types(); ++t) {
            for (idx_t p : build_halo.periodic_cells_local_index_[t]) {
                cells_to_edit[i++] = c + p;
            }
            c += cells.elements(t).size();
        }
        glb_idx_max = cells.global_index().metadata().getLong("max", 0);
    }

    std::vector<gidx_t> glb_idx(cells_to_edit.size());
    const int nb_cells = static_cast<int>(glb_idx.size());
    for (int i = 0; i < nb_cells; ++i) {
        glb_idx[i] = cells_glb_idx(cells_to_edit[i]);
    }

    // 1) Gather all global indices, together with location

    std::vector<int> recvcounts(mpi::size());
    std::vector<int> recvdispls(mpi::size());

    ATLAS_TRACE_MPI(GATHER) { mpi::comm().gather(nb_cells, recvcounts, root); }
    idx_t glb_nb_cells = std::accumulate(recvcounts.begin(), recvcounts.end(), 0);

    recvdispls[0] = 0;
    for (int jpart = 1; jpart < nparts; ++jpart) {  // start at 1
        recvdispls[jpart] = recvcounts[jpart - 1] + recvdispls[jpart - 1];
    }

    std::vector<gidx_t> glb_idx_gathered(glb_nb_cells);
    ATLAS_TRACE_MPI(GATHER) {
        mpi::comm().gatherv(glb_idx.data(), glb_idx.size(), glb_idx_gathered.data(), recvcounts.data(),
                            recvdispls.data(), root);
    }

    // 2) Sort all global indices, and renumber from 1 to glb_nb_edges
    std::vector<Entity> cell_sort;
    cell_sort.reserve(glb_nb_cells);
    for (idx_t jcell = 0; jcell < glb_nb_cells; ++jcell) {
        cell_sort.emplace_back(glb_idx_gathered[jcell], jcell);
    }

    ATLAS_TRACE_SCOPE("sort on rank 0") { std::sort(cell_sort.begin(), cell_sort.end()); }

    gidx_t gid = glb_idx_max + 1;
    for (idx_t jcell = 0; jcell < glb_nb_cells; ++jcell) {
        if (jcell > 0 && cell_sort[jcell].g != cell_sort[jcell - 1].g) {
            ++gid;
        }
        idx_t icell             = cell_sort[jcell].i;
        glb_idx_gathered[icell] = gid;
    }

    // 3) Scatter renumbered back
    ATLAS_TRACE_MPI(SCATTER) {
        mpi::comm().scatterv(glb_idx_gathered.data(), recvcounts.data(), recvdispls.data(), glb_idx.data(),
                             static_cast<int>(glb_idx.size()), root);
    }

    for (int jcell = 0; jcell < nb_cells; ++jcell) {
        cells_glb_idx(cells_to_edit[jcell]) = glb_idx[jcell];
    }

    // nodes_glb_idx.dump( Log::info() );
    // Log::info() << std::endl;
    cells.global_index().metadata().set("human_readable", true);
}

// -------------------------------------------------------------------------------------

using uid_t = gidx_t;

// ------------------------------------------------------------------
class BuildHaloHelper;

void increase_halo(Mesh& mesh);
void increase_halo_interior(BuildHaloHelper&);

class EastWest : public util::PeriodicTransform {
public:
    EastWest() { x_translation_ = -360.; }
};

class WestEast : public util::PeriodicTransform {
public:
    WestEast() { x_translation_ = 360.; }
};

using Node2Elem = std::vector<std::vector<idx_t>>;

void build_lookup_node2elem(const Mesh& mesh, Node2Elem& node2elem) {
    ATLAS_TRACE();

    const mesh::Nodes& nodes = mesh.nodes();

    node2elem.resize(nodes.size());
    for (idx_t jnode = 0; jnode < nodes.size(); ++jnode) {
        node2elem[jnode].clear();
        node2elem[jnode].reserve(12);
    }

    const mesh::HybridElements::Connectivity& elem_nodes = mesh.cells().node_connectivity();
    auto field_flags                                     = array::make_view<int, 1>(mesh.cells().flags());
    auto patched                                         = [&field_flags](idx_t e) {
        using Topology = atlas::mesh::Nodes::Topology;
        return Topology::check(field_flags(e), Topology::PATCH);
    };

    idx_t nb_elems = mesh.cells().size();
    for (idx_t elem = 0; elem < nb_elems; ++elem) {
        if (not patched(elem)) {
            for (idx_t n = 0; n < elem_nodes.cols(elem); ++n) {
                int node = elem_nodes(elem, n);
                node2elem[node].push_back(elem);
            }
        }
    }
}

void accumulate_partition_bdry_nodes_old(Mesh& mesh, std::vector<idx_t>& bdry_nodes) {
    ATLAS_TRACE();

    std::set<idx_t> bdry_nodes_set;

    std::vector<idx_t> facet_nodes;
    std::vector<idx_t> connectivity_facet_to_elem;

    facet_nodes.reserve(mesh.nodes().size() * 4);
    connectivity_facet_to_elem.reserve(facet_nodes.capacity() * 2);

    idx_t nb_facets(0);
    idx_t nb_inner_facets(0);
    idx_t missing_value;
    accumulate_facets(
        /*in*/ mesh.cells(),
        /*in*/ mesh.nodes(),
        /*out*/ facet_nodes,  // shape(nb_facets,nb_nodes_per_facet)
        /*out*/ connectivity_facet_to_elem,
        /*out*/ nb_facets,
        /*out*/ nb_inner_facets,
        /*out*/ missing_value);

    for (idx_t jface = 0; jface < nb_facets; ++jface) {
        if (connectivity_facet_to_elem[jface * 2 + 1] == missing_value) {
            for (idx_t jnode = 0; jnode < 2; ++jnode)  // 2 nodes per face
            {
                bdry_nodes_set.insert(facet_nodes[jface * 2 + jnode]);
            }
        }
    }
    bdry_nodes = std::vector<idx_t>(bdry_nodes_set.begin(), bdry_nodes_set.end());
}

void accumulate_partition_bdry_nodes(Mesh& mesh, idx_t halo, std::vector<idx_t>& bdry_nodes) {
#ifndef ATLAS_103
    /* deprecated */
    accumulate_partition_bdry_nodes_old(mesh, bdry_nodes);
#else
    ATLAS_TRACE();
    const Mesh::Polygon& polygon = mesh.polygon(halo);
    bdry_nodes                   = std::vector<int>(polygon.begin(), polygon.end());
#endif

#ifdef ATLAS_103_SORT
    /* not required */
    std::sort(bdry_nodes.begin(), bdry_nodes.end());
#endif
}

template <typename Predicate>
std::vector<idx_t> filter_nodes(std::vector<idx_t> nodes, const Predicate& predicate) {
    std::vector<idx_t> filtered;
    filtered.reserve(nodes.size());
    for (idx_t inode : nodes) {
        if (predicate(inode)) {
            filtered.push_back(inode);
        }
    }
    return filtered;
}

class Notification {
public:
    void add_error(const std::string& note, const eckit::CodeLocation& loc) {
        notes.push_back(note + " @ " + std::string(loc));
    }
    void add_error(const std::string& note) { notes.push_back(note); }

    bool error() const { return notes.size() > 0; }
    void reset() { notes.clear(); }

    std::string str() const {
        std::stringstream stream;
        for (size_t jnote = 0; jnote < notes.size(); ++jnote) {
            if (jnote > 0) {
                stream << "\n";
            }
            stream << notes[jnote];
        }
        return stream.str();
    }

    operator std::string() const { return str(); }

private:
    friend std::ostream& operator<<(std::ostream& s, const Notification& notes) {
        s << notes.str();
        return s;
    }

private:
    std::vector<std::string> notes;
};

using Uid2Node = std::unordered_map<uid_t, idx_t>;


void build_lookup_uid2node(Mesh& mesh, Uid2Node& uid2node) {
    ATLAS_TRACE();
    Notification notes;
    mesh::Nodes& nodes                  = mesh.nodes();
    array::ArrayView<double, 2> xy      = array::make_view<double, 2>(nodes.xy());
    array::ArrayView<gidx_t, 1> glb_idx = array::make_view<gidx_t, 1>(nodes.global_index());
    idx_t nb_nodes                      = nodes.size();

    UniqueLonLat compute_uid(mesh);

    uid2node.clear();
    for (idx_t jnode = 0; jnode < nb_nodes; ++jnode) {
        uid_t uid     = compute_uid(jnode);
        bool inserted = uid2node.insert(std::make_pair(uid, jnode)).second;
        if (not inserted) {
                int other = uid2node[uid];
                std::stringstream msg;
                msg << std::setprecision(10) << std::fixed << "Node uid: " << uid << "   " << glb_idx(jnode) << " xy("
                    << xy(jnode, XX) << "," << xy(jnode, YY) << ")";
                if (nodes.has_field("ij")) {
                    auto ij = array::make_view<idx_t, 2>(nodes.field("ij"));
                    msg << " ij(" << ij(jnode, XX) << "," << ij(jnode, YY) << ")";
                }
                msg << " has already been added as node " << glb_idx(other) << " (" << xy(other, XX) << "," << xy(other, YY)
                    << ")";
                if (nodes.has_field("ij")) {
                    auto ij = array::make_view<idx_t, 2>(nodes.field("ij"));
                    msg << " ij(" << ij(other, XX) << "," << ij(other, YY) << ")";
                }
                notes.add_error(msg.str());
        }
    }
    if (notes.error()) {
        throw_Exception(notes.str(), Here());
    }
}

/// For given nodes, find new connected elements
void accumulate_elements(const Mesh& mesh, const mpi::BufferView<uid_t>& request_node_uid, const Uid2Node& uid2node,
                         const Node2Elem& node2elem, std::vector<idx_t>& found_elements,
                         std::set<uid_t>& new_nodes_uid) {

    ATLAS_TRACE();
    const mesh::HybridElements::Connectivity& elem_nodes = mesh.cells().node_connectivity();
    const auto elem_part                                 = array::make_view<int, 1>(mesh.cells().partition());

    const idx_t nb_nodes         = mesh.nodes().size();
    const idx_t nb_request_nodes = static_cast<idx_t>(request_node_uid.size());
    const int mpi_rank           = static_cast<int>(mpi::rank());

    std::unordered_set<idx_t> found_elements_set;
    found_elements_set.reserve(nb_request_nodes*2);

    for (idx_t jnode = 0; jnode < nb_request_nodes; ++jnode) {
        uid_t uid = request_node_uid(jnode);

        idx_t inode = -1;
        // search and get node index for uid
        auto found = uid2node.find(uid);
        if (found != uid2node.end()) {
            inode = found->second;
        }
        if (inode != -1 && inode < nb_nodes) {
            for (const idx_t e : node2elem[inode]) {
                if (elem_part(e) == mpi_rank) {
                    found_elements_set.insert(e);
                }
            }
        }
    }

    // found_bdry_elements_set now contains elements for the nodes
    found_elements = std::vector<idx_t>(found_elements_set.begin(), found_elements_set.end());

    UniqueLonLat compute_uid(mesh);

    // Collect all nodes
    new_nodes_uid.clear();
    for (const idx_t e : found_elements) {
        idx_t nb_elem_nodes = elem_nodes.cols(e);
        for (idx_t n = 0; n < nb_elem_nodes; ++n) {
            new_nodes_uid.insert(compute_uid(elem_nodes(e, n)));
        }
    }

    // Remove nodes we already have in the request-buffer
    for (idx_t jnode = 0; jnode < nb_request_nodes; ++jnode) {
        new_nodes_uid.erase(request_node_uid(jnode));
    }
}

class BuildHaloHelper {
public:
    struct Buffers {
        std::vector<std::vector<int>> node_part;

        std::vector<std::vector<int>> node_ridx;

        std::vector<std::vector<int>> node_flags;

        std::vector<std::vector<uid_t>> node_glb_idx;

        std::vector<std::vector<double>> node_xy;

        std::vector<std::vector<uid_t>> elem_glb_idx;

        std::vector<std::vector<uid_t>> elem_nodes_id;

        std::vector<std::vector<int>> elem_nodes_displs;

        std::vector<std::vector<int>> elem_part;

        std::vector<std::vector<int>> elem_ridx;

        std::vector<std::vector<int>> elem_flags;

        std::vector<std::vector<int>> elem_type;

        Buffers(Mesh&) {
            const idx_t mpi_size = static_cast<idx_t>(mpi::size());

            node_part.resize(mpi_size);
            node_ridx.resize(mpi_size);
            node_flags.resize(mpi_size);
            node_glb_idx.resize(mpi_size);
            node_xy.resize(mpi_size);
            elem_glb_idx.resize(mpi_size);
            elem_nodes_id.resize(mpi_size);
            elem_nodes_displs.resize(mpi_size);
            elem_part.resize(mpi_size);
            elem_ridx.resize(mpi_size);
            elem_flags.resize(mpi_size);
            elem_type.resize(mpi_size);
        }

        void print(std::ostream& os) const {
            const idx_t mpi_size = static_cast<idx_t>(mpi::size());
            os << "Nodes\n"
               << "-----\n";
            idx_t n(0);
            for (idx_t jpart = 0; jpart < mpi_size; ++jpart) {
                for (idx_t jnode = 0; jnode < node_glb_idx[jpart].size(); ++jnode ) {
                    auto g = node_glb_idx[jpart][jnode];
                    auto p = node_part[jpart][jnode];
                    os << std::setw(4) << n++ << " : [ p" << p << " ]  " << g << "\n";
                }
            }
            os << std::flush;
            os << "Cells\n"
               << "-----\n";
            idx_t e(0);
            for (idx_t jpart = 0; jpart < mpi_size; ++jpart) {
                const idx_t nb_elem = static_cast<idx_t>(elem_glb_idx[jpart].size());
                for (idx_t jelem = 0; jelem < nb_elem; ++jelem) {
                    os << std::setw(4) << e++ << " :  [ t" << elem_type[jpart][jelem] << " -- p"
                       << elem_part[jpart][jelem] << "]  " << elem_glb_idx[jpart][jelem] << "\n";
                }
            }
            os << std::flush;
        }
        friend std::ostream& operator<<(std::ostream& out, const Buffers& b) {
            b.print(out);
            return out;
        }
    };

    static void all_to_all(Buffers& send, Buffers& recv) {
        ATLAS_TRACE();
        const eckit::mpi::Comm& comm = mpi::comm();

        ATLAS_TRACE_MPI(ALLTOALL) {
            comm.allToAll(send.node_glb_idx, recv.node_glb_idx);
            comm.allToAll(send.node_part, recv.node_part);
            comm.allToAll(send.node_ridx, recv.node_ridx);
            comm.allToAll(send.node_flags, recv.node_flags);
            comm.allToAll(send.node_xy, recv.node_xy);
            comm.allToAll(send.elem_glb_idx, recv.elem_glb_idx);
            comm.allToAll(send.elem_nodes_id, recv.elem_nodes_id);
            comm.allToAll(send.elem_part, recv.elem_part);
            comm.allToAll(send.elem_ridx, recv.elem_ridx);
            comm.allToAll(send.elem_type, recv.elem_type);
            comm.allToAll(send.elem_flags, recv.elem_flags);
            comm.allToAll(send.elem_nodes_displs, recv.elem_nodes_displs);
        }
    }

    struct Status {
        std::vector<idx_t> new_periodic_ghost_points;
        std::vector<std::vector<idx_t>> new_periodic_ghost_cells;
    } status;

public:
    BuildHalo& builder_;
    Mesh& mesh;
    array::ArrayView<double, 2> xy;
    array::ArrayView<double, 2> lonlat;
    array::ArrayView<gidx_t, 1> glb_idx;
    array::ArrayView<int, 1> part;
    array::IndexView<idx_t, 1> ridx;
    array::ArrayView<int, 1> flags;
    array::ArrayView<int, 1> halo;
    array::ArrayView<int, 1> ghost;
    mesh::HybridElements::Connectivity* elem_nodes;
    array::ArrayView<int, 1> elem_part;
    array::IndexView<idx_t, 1> elem_ridx;
    array::ArrayView<int, 1> elem_flags;
    array::ArrayView<gidx_t, 1> elem_glb_idx;

    std::vector<idx_t> bdry_nodes;
    Node2Elem node_to_elem;
    Uid2Node uid2node;
    UniqueLonLat compute_uid;
    idx_t halosize;

public:
    BuildHaloHelper(BuildHalo& builder, Mesh& _mesh):
        builder_(builder),
        mesh(_mesh),
        xy(array::make_view<double, 2>(mesh.nodes().xy())),
        lonlat(array::make_view<double, 2>(mesh.nodes().lonlat())),
        glb_idx(array::make_view<gidx_t, 1>(mesh.nodes().global_index())),
        part(array::make_view<int, 1>(mesh.nodes().partition())),
        ridx(array::make_indexview<idx_t, 1>(mesh.nodes().remote_index())),
        flags(array::make_view<int, 1>(mesh.nodes().flags())),
        halo(array::make_view<int, 1>(mesh.nodes().halo())),
        ghost(array::make_view<int, 1>(mesh.nodes().ghost())),
        elem_nodes(&mesh.cells().node_connectivity()),
        elem_part(array::make_view<int, 1>(mesh.cells().partition())),
        elem_ridx(array::make_indexview<idx_t, 1>(mesh.cells().remote_index())),
        elem_flags(array::make_view<int, 1>(mesh.cells().flags())),
        elem_glb_idx(array::make_view<gidx_t, 1>(mesh.cells().global_index())),
        compute_uid(mesh) {
        halosize = 0;
        mesh.metadata().get("halo", halosize);
        // update();
    }

    void update() {
        compute_uid.update();
        mesh::Nodes& nodes = mesh.nodes();
        xy                 = array::make_view<double, 2>(nodes.xy());
        lonlat             = array::make_view<double, 2>(nodes.lonlat());
        glb_idx            = array::make_view<gidx_t, 1>(nodes.global_index());
        part               = array::make_view<int, 1>(nodes.partition());
        ridx               = array::make_indexview<idx_t, 1>(nodes.remote_index());
        flags              = array::make_view<int, 1>(nodes.flags());
        ghost              = array::make_view<int, 1>(nodes.ghost());
        halo               = array::make_view<int, 1>(nodes.halo());

        elem_nodes   = &mesh.cells().node_connectivity();
        elem_part    = array::make_view<int, 1>(mesh.cells().partition());
        elem_ridx    = array::make_indexview<idx_t, 1>(mesh.cells().remote_index());
        elem_flags   = array::make_view<int, 1>(mesh.cells().flags());
        elem_glb_idx = array::make_view<gidx_t, 1>(mesh.cells().global_index());
    }

    template <typename NodeContainer, typename ElementContainer>
    void fill_sendbuffer(Buffers& buf, const NodeContainer& nodes_uid, const ElementContainer& elems, const int p) {
        ATLAS_TRACE();

        idx_t nb_nodes = static_cast<idx_t>(nodes_uid.size());
        buf.node_glb_idx[p].resize(nb_nodes);
        buf.node_part[p].resize(nb_nodes);
        buf.node_ridx[p].resize(nb_nodes);
        buf.node_flags[p].resize(nb_nodes, Topology::NONE);
        buf.node_xy[p].resize(2 * nb_nodes);

        idx_t jnode = 0;
        typename NodeContainer::const_iterator it;
        for (it = nodes_uid.begin(); it != nodes_uid.end(); ++it, ++jnode) {
            uid_t uid = *it;
            auto found = uid2node.find(uid);
            if (found != uid2node.end())  // Point exists inside domain
            {
                idx_t node                     = found->second;
                buf.node_glb_idx[p][jnode]     = glb_idx(node);
                buf.node_part[p][jnode]        = part(node);
                buf.node_ridx[p][jnode]        = ridx(node);
                buf.node_xy[p][jnode * 2 + XX] = xy(node, XX);
                buf.node_xy[p][jnode * 2 + YY] = xy(node, YY);
                Topology::set(buf.node_flags[p][jnode], flags(node) | Topology::GHOST);
            }
            else {
                Log::warning() << "Node with uid " << uid << " needed by [" << p << "] was not found in ["
                               << mpi::rank() << "]." << std::endl;
                ATLAS_ASSERT(false);
            }
        }

        idx_t nb_elems = static_cast<idx_t>(elems.size());

        idx_t nb_elem_nodes(0);
        for (idx_t jelem = 0; jelem < nb_elems; ++jelem) {
            idx_t ielem = elems[jelem];
            nb_elem_nodes += elem_nodes->cols(ielem);
        }

        buf.elem_glb_idx[p].resize(nb_elems);
        buf.elem_part[p].resize(nb_elems);
        buf.elem_ridx[p].resize(nb_elems);
        buf.elem_flags[p].resize(nb_elems, Topology::NONE);
        buf.elem_type[p].resize(nb_elems);
        buf.elem_nodes_id[p].resize(nb_elem_nodes);
        buf.elem_nodes_displs[p].resize(nb_elems);
        idx_t jelemnode(0);
        for (idx_t jelem = 0; jelem < nb_elems; ++jelem) {
            buf.elem_nodes_displs[p][jelem] = jelemnode;
            idx_t ielem                     = elems[jelem];

            buf.elem_glb_idx[p][jelem] = elem_glb_idx(ielem);
            buf.elem_part[p][jelem]    = elem_part(ielem);
            buf.elem_ridx[p][jelem]    = elem_ridx(ielem);
            Topology::set(buf.elem_flags[p][jelem], elem_flags(ielem));
            buf.elem_type[p][jelem] = mesh.cells().type_idx(ielem);
            for (idx_t jnode = 0; jnode < elem_nodes->cols(ielem); ++jnode) {
                buf.elem_nodes_id[p][jelemnode++] = compute_uid((*elem_nodes)(ielem, jnode));
            }
            Topology::set(buf.elem_flags[p][jelem], Topology::GHOST);
        }
    }

    template <typename NodeContainer, typename ElementContainer>
    void fill_sendbuffer(Buffers& buf, const NodeContainer& nodes_uid, const ElementContainer& elems,
                         const util::PeriodicTransform& transform, int newflags, const int p) {
        // ATLAS_TRACE();

        idx_t nb_nodes = static_cast<idx_t>(nodes_uid.size());
        buf.node_glb_idx[p].resize(nb_nodes);
        buf.node_part[p].resize(nb_nodes);
        buf.node_ridx[p].resize(nb_nodes);
        buf.node_flags[p].resize(nb_nodes, Topology::NONE);
        buf.node_xy[p].resize(2 * nb_nodes);

        int jnode = 0;
        typename NodeContainer::const_iterator it;
        for (it = nodes_uid.begin(); it != nodes_uid.end(); ++it, ++jnode) {
            uid_t uid = *it;

            auto found = uid2node.find(uid);
            if (found != uid2node.end())  // Point exists inside domain
            {
                int node                       = found->second;
                buf.node_part[p][jnode]        = part(node);
                buf.node_ridx[p][jnode]        = ridx(node);
                buf.node_xy[p][jnode * 2 + XX] = xy(node, XX);
                buf.node_xy[p][jnode * 2 + YY] = xy(node, YY);
                transform(&buf.node_xy[p][jnode * 2], -1);
                // Global index of node is based on UID of destination
                buf.node_glb_idx[p][jnode] = util::unique_lonlat(&buf.node_xy[p][jnode * 2]);
                Topology::set(buf.node_flags[p][jnode], newflags);
            }
            else {
                Log::warning() << "Node with uid " << uid << " needed by [" << p << "] was not found in ["
                               << mpi::rank() << "]." << std::endl;
                ATLAS_ASSERT(false);
            }
        }

        idx_t nb_elems = static_cast<idx_t>(elems.size());

        idx_t nb_elem_nodes(0);
        for (idx_t jelem = 0; jelem < nb_elems; ++jelem) {
            idx_t ielem = elems[jelem];
            nb_elem_nodes += elem_nodes->cols(ielem);
        }

        buf.elem_glb_idx[p].resize(nb_elems);
        buf.elem_part[p].resize(nb_elems);
        buf.elem_ridx[p].resize(nb_elems);
        buf.elem_flags[p].resize(nb_elems, Topology::NONE);
        buf.elem_type[p].resize(nb_elems);
        buf.elem_nodes_id[p].resize(nb_elem_nodes);
        buf.elem_nodes_displs[p].resize(nb_elems);
        idx_t jelemnode(0);
        for (idx_t jelem = 0; jelem < nb_elems; ++jelem) {
            buf.elem_nodes_displs[p][jelem] = jelemnode;
            idx_t ielem                     = elems[jelem];
            buf.elem_part[p][jelem]         = elem_part(ielem);
            buf.elem_ridx[p][jelem]         = elem_ridx(ielem);
            Topology::set(buf.elem_flags[p][jelem], elem_flags(ielem) | newflags);
            buf.elem_type[p][jelem] = mesh.cells().type_idx(ielem);
            std::vector<double> crds(elem_nodes->cols(ielem) * 2);
            for (idx_t jnode = 0; jnode < elem_nodes->cols(ielem); ++jnode) {
                double crd[] = {xy((*elem_nodes)(ielem, jnode), XX), xy((*elem_nodes)(ielem, jnode), YY)};
                transform(crd, -1);
                buf.elem_nodes_id[p][jelemnode++] = util::unique_lonlat(crd);
                crds[jnode * 2 + XX]              = crd[XX];
                crds[jnode * 2 + YY]              = crd[YY];
            }
            // Global index of element is based on UID of destination

            buf.elem_glb_idx[p][jelem] = -util::unique_lonlat(crds.data(), elem_nodes->cols(ielem));
        }
    }

    void add_nodes(Buffers& buf) {
        ATLAS_TRACE();

        const idx_t mpi_size = static_cast<idx_t>(mpi::size());

        mesh::Nodes& nodes = mesh.nodes();
        int nb_nodes       = nodes.size();

        // Nodes might be duplicated from different Tasks. We need to identify
        // unique entries
        std::vector<uid_t> node_uid(nb_nodes);
        std::set<uid_t> new_node_uid;
        {
            ATLAS_TRACE("compute node_uid");
            for (int jnode = 0; jnode < nb_nodes; ++jnode) {
                node_uid[jnode] = compute_uid(jnode);
            }
            std::sort(node_uid.begin(), node_uid.end());
        }
        auto node_already_exists = [&node_uid, &new_node_uid](uid_t uid) {
            std::vector<uid_t>::iterator it = std::lower_bound(node_uid.begin(), node_uid.end(), uid);
            bool not_found                  = (it == node_uid.end() || uid < *it);
            if (not_found) {
                bool inserted = new_node_uid.insert(uid).second;
                return not inserted;
            }
            else {
                return true;
            }
        };

        std::vector<std::vector<int>> rfn_idx(mpi_size);
        for (idx_t jpart = 0; jpart < mpi_size; ++jpart) {
            rfn_idx[jpart].reserve(buf.node_glb_idx[jpart].size());
        }

        int nb_new_nodes = 0;
        for (idx_t jpart = 0; jpart < mpi_size; ++jpart) {
            const idx_t nb_nodes_on_part = static_cast<idx_t>(buf.node_glb_idx[jpart].size());
            for (idx_t n = 0; n < nb_nodes_on_part; ++n) {
                std::array<double,2> crd{buf.node_xy[jpart][n * 2 + XX], buf.node_xy[jpart][n * 2 + YY]};
                mesh.projection().xy2lonlat(crd.data());

                if (not node_already_exists(util::unique_lonlat(crd))) {
                    rfn_idx[jpart].push_back(n);
                }
            }
            nb_new_nodes += rfn_idx[jpart].size();
        }

        // Resize nodes
        // ------------
        nodes.resize(nb_nodes + nb_new_nodes);
        flags   = array::make_view<int, 1>(nodes.flags());
        halo    = array::make_view<int, 1>(nodes.halo());
        glb_idx = array::make_view<gidx_t, 1>(nodes.global_index());
        part    = array::make_view<int, 1>(nodes.partition());
        ridx    = array::make_indexview<idx_t, 1>(nodes.remote_index());
        xy      = array::make_view<double, 2>(nodes.xy());
        lonlat  = array::make_view<double, 2>(nodes.lonlat());
        ghost   = array::make_view<int, 1>(nodes.ghost());

        compute_uid.update();

        // Add new nodes
        // -------------
        int new_node = 0;
        for (idx_t jpart = 0; jpart < mpi_size; ++jpart) {
            for (size_t n = 0; n < rfn_idx[jpart].size(); ++n) {
                int loc_idx   = nb_nodes + new_node;
                halo(loc_idx) = halosize + 1;
                Topology::reset(flags(loc_idx), buf.node_flags[jpart][rfn_idx[jpart][n]]);
                ghost(loc_idx)   = Topology::check(flags(loc_idx), Topology::GHOST);
                glb_idx(loc_idx) = buf.node_glb_idx[jpart][rfn_idx[jpart][n]];
                part(loc_idx)    = buf.node_part[jpart][rfn_idx[jpart][n]];
                ridx(loc_idx)    = buf.node_ridx[jpart][rfn_idx[jpart][n]];
                PointXY pxy(&buf.node_xy[jpart][rfn_idx[jpart][n] * 2]);
                xy(loc_idx, XX)     = pxy.x();
                xy(loc_idx, YY)     = pxy.y();
                PointLonLat pll     = mesh.projection().lonlat(pxy);
                lonlat(loc_idx, XX) = pll.lon();
                lonlat(loc_idx, YY) = pll.lat();

                if (Topology::check(flags(loc_idx), Topology::PERIODIC) and
                    not Topology::check(flags(loc_idx), Topology::BC)) {
                    status.new_periodic_ghost_points.push_back(loc_idx);
                }

                // make sure new node was not already there
                {
                    uid_t uid                = compute_uid(loc_idx);
                    auto found = uid2node.find(uid);
                    if (found != uid2node.end()) {
                        int other = found->second;
                        std::stringstream msg;
                        msg << "New node loc " << loc_idx << " with uid " << uid << ":\n"
                            << glb_idx(loc_idx) << "(" << xy(loc_idx, XX) << "," << xy(loc_idx, YY) << ")\n";
                        msg << "Existing already loc " << other << "  :  " << glb_idx(other) << "(" << xy(other, XX)
                            << "," << xy(other, YY) << ")\n";
                        throw_Exception(msg.str(), Here());
                    }
                    uid2node[uid] = nb_nodes + new_node;
                }
                ++new_node;
            }
        }
    }

    void add_elements(Buffers& buf) {
        ATLAS_TRACE();

        const idx_t mpi_size = static_cast<idx_t>(mpi::size());
        auto cell_gidx       = array::make_view<gidx_t, 1>(mesh.cells().global_index());
        // Elements might be duplicated from different Tasks. We need to identify
        // unique entries
        int nb_elems = mesh.cells().size();
        //    std::set<uid_t> elem_uid;
        std::vector<uid_t> elem_uid(2 * nb_elems);
        std::set<uid_t> new_elem_uid;
        {
            ATLAS_TRACE("compute elem_uid");
            atlas_omp_parallel_for (int jelem = 0; jelem < nb_elems; ++jelem) {
                elem_uid[jelem * 2 + 0] = -compute_uid(elem_nodes->row(jelem));
                elem_uid[jelem * 2 + 1] = cell_gidx(jelem);
            }
            omp::sort(elem_uid.begin(), elem_uid.end());
        }
        auto element_already_exists = [&elem_uid, &new_elem_uid](uid_t uid) -> bool {
            std::vector<uid_t>::iterator it = std::lower_bound(elem_uid.begin(), elem_uid.end(), uid);
            bool not_found                  = (it == elem_uid.end() || uid < *it);
            if (not_found) {
                bool inserted = new_elem_uid.insert(uid).second;
                return not inserted;
            }
            else {
                return true;
            }
        };

        if (not status.new_periodic_ghost_cells.size()) {
            status.new_periodic_ghost_cells.resize(mesh.cells().nb_types());
        }

        std::vector<std::vector<idx_t>> received_new_elems(mpi_size);
        for (idx_t jpart = 0; jpart < mpi_size; ++jpart) {
            received_new_elems[jpart].reserve(buf.elem_glb_idx[jpart].size());
        }

        idx_t nb_new_elems(0);
        for (idx_t jpart = 0; jpart < mpi_size; ++jpart) {
            const idx_t nb_elems_from_part = static_cast<idx_t>(buf.elem_glb_idx[jpart].size());
            for (idx_t e = 0; e < nb_elems_from_part; ++e) {
                if (element_already_exists(buf.elem_glb_idx[jpart][e]) == false) {
                    received_new_elems[jpart].emplace_back(e);
                }
            }
            nb_new_elems += received_new_elems[jpart].size();
        }

        std::vector<std::vector<std::vector<int>>> elements_of_type(mesh.cells().nb_types(),
                                                                    std::vector<std::vector<int>>(mpi_size));
        std::vector<idx_t> nb_elements_of_type(mesh.cells().nb_types(), 0);

        for (idx_t jpart = 0; jpart < mpi_size; ++jpart) {
            for (const idx_t ielem : received_new_elems[jpart]) {
                elements_of_type[buf.elem_type[jpart][ielem]][jpart].emplace_back(ielem);
                ++nb_elements_of_type[buf.elem_type[jpart][ielem]];
            }
        }

        for (idx_t t = 0; t < mesh.cells().nb_types(); ++t) {
            const std::vector<std::vector<int>>& elems = elements_of_type[t];
            mesh::Elements& elements                   = mesh.cells().elements(t);

            // Add new elements
            BlockConnectivity& node_connectivity = elements.node_connectivity();
            if (nb_elements_of_type[t] == 0) {
                continue;
            }

            idx_t old_size      = elements.size();
            idx_t new_elems_pos = elements.add(nb_elements_of_type[t]);

            auto elem_type_glb_idx = elements.view<gidx_t, 1>(mesh.cells().global_index());
            auto elem_type_part    = elements.view<int, 1>(mesh.cells().partition());
            auto elem_type_ridx    = elements.indexview<idx_t, 1>(mesh.cells().remote_index());
            auto elem_type_halo    = elements.view<int, 1>(mesh.cells().halo());
            auto elem_type_flags   = elements.view<int, 1>(mesh.cells().flags());

            // Copy information in new elements
            idx_t new_elem(0);
            for (idx_t jpart = 0; jpart < mpi_size; ++jpart) {
                for (const idx_t jelem : elems[jpart]) {
                    int loc_idx                = new_elems_pos + new_elem;
                    elem_type_glb_idx(loc_idx) = std::abs(buf.elem_glb_idx[jpart][jelem]);
                    elem_type_part(loc_idx)    = buf.elem_part[jpart][jelem];
                    elem_type_ridx(loc_idx)    = buf.elem_ridx[jpart][jelem];
                    elem_type_halo(loc_idx)    = halosize + 1;
                    elem_type_flags(loc_idx)   = buf.elem_flags[jpart][jelem];
                    for (idx_t n = 0; n < node_connectivity.cols(); ++n) {
                        node_connectivity.set(
                            loc_idx, n, uid2node[buf.elem_nodes_id[jpart][buf.elem_nodes_displs[jpart][jelem] + n]]);
                    }

                    if (Topology::check(elem_type_flags(loc_idx), Topology::PERIODIC)) {
                        status.new_periodic_ghost_cells[t].push_back(old_size + new_elem);
                    }
                    ++new_elem;
                }
            }
        }
    }

    std::vector<std::vector<idx_t>> gather_nb_elements() {
        idx_t nb_types = mesh.cells().nb_types();
        idx_t mpi_size = mpi::comm().size();
        std::vector<idx_t> elems_per_type(nb_types);
        for (idx_t t = 0; t < nb_types; ++t) {
            elems_per_type[t] = mesh.cells().elements(t).size();
        }
        atlas::mpi::Buffer<idx_t, 1> recv(mpi_size);
        ATLAS_TRACE_MPI(ALLGATHER) { mpi::comm().allGatherv(elems_per_type.begin(), elems_per_type.end(), recv); }

        std::vector<std::vector<idx_t>> result(nb_types, std::vector<idx_t>(mpi_size));
        for (idx_t p = 0; p < mpi_size; ++p) {
            auto recv_p = recv[p];
            for (idx_t t = 0; t < nb_types; ++t) {
                result[t][p] = recv_p[t];
            }
        }
        return result;
    };


    void add_buffers(Buffers& buf) {
        add_nodes(buf);

        auto nb_elements_pre = gather_nb_elements();

        add_elements(buf);

        auto nb_elements_post = gather_nb_elements();

        // Renumber remote_index as the number of elements per type might have grown
        {
            auto nb_partitions = mpi::comm().size();
            auto nb_types      = mesh.cells().nb_types();

            std::vector<std::vector<idx_t>> accumulated_diff(nb_types, std::vector<idx_t>(nb_partitions));
            for (idx_t t = 1; t < nb_types; ++t) {
                for (idx_t p = 0; p < nb_partitions; ++p) {
                    accumulated_diff[t][p] =
                        accumulated_diff[t - 1][p] + nb_elements_post[t - 1][p] - nb_elements_pre[t - 1][p];
                }
            }

            for (idx_t t = 0; t < nb_types; ++t) {
                auto& elements = mesh.cells().elements(t);
                auto partition = elements.view<int, 1>(mesh.cells().partition());
                auto ridx      = elements.indexview<idx_t, 1>(mesh.cells().remote_index());

                auto& accumulated_diff_t = accumulated_diff[t];

                for (idx_t elem = 0; elem < elements.size(); ++elem) {
                    ridx(elem) += accumulated_diff_t[partition(elem)];
                }
            }
        }

        update();
    }
};

namespace {
void gather_bdry_nodes(const BuildHaloHelper& helper, const std::vector<uid_t>& send,
                       atlas::mpi::Buffer<uid_t, 1>& recv, bool periodic = false) {
    auto& comm = mpi::comm();
#ifndef ATLAS_103
    /* deprecated */
    ATLAS_TRACE("gather_bdry_nodes old way");
    {
        ATLAS_TRACE_MPI(ALLGATHER) { comm.allGatherv(send.begin(), send.end(), recv); }
    }
#else
    ATLAS_TRACE();

    Mesh::PartitionGraph::Neighbours neighbours = helper.mesh.nearestNeighbourPartitions();
    if (periodic) {
        // add own rank to neighbours to allow periodicity with self (pole caps)
        idx_t rank = comm.rank();
        neighbours.insert(std::upper_bound(neighbours.begin(), neighbours.end(), rank), rank);
    }

    const idx_t mpi_size = comm.size();
    const int counts_tag = 0;
    const int buffer_tag = 1;

    std::vector<eckit::mpi::Request> counts_requests;
    counts_requests.reserve(neighbours.size());
    std::vector<eckit::mpi::Request> buffer_requests;
    buffer_requests.reserve(neighbours.size());

    int sendcnt = send.size();
    ATLAS_TRACE_MPI(ISEND) {
        for (idx_t to : neighbours) {
            counts_requests.push_back(comm.iSend(sendcnt, to, counts_tag));
        }
    }

    recv.counts.assign(0, mpi_size);

    ATLAS_TRACE_MPI(IRECEIVE) {
        for (idx_t from : neighbours) {
            counts_requests.push_back(comm.iReceive(recv.counts[from], from, counts_tag));
        }
    }

    ATLAS_TRACE_MPI(ISEND) {
        for (idx_t to : neighbours) {
            buffer_requests.push_back(comm.iSend(send.data(), send.size(), to, buffer_tag));
        }
    }

    ATLAS_TRACE_MPI(WAIT) {
        for (auto request : counts_requests) {
            comm.wait(request);
        }
    }

    recv.displs[0] = 0;
    recv.cnt       = recv.counts[0];
    for (idx_t jpart = 1; jpart < mpi_size; ++jpart) {
        recv.displs[jpart] = recv.displs[jpart - 1] + recv.counts[jpart - 1];
        recv.cnt += recv.counts[jpart];
    }
    recv.buffer.resize(recv.cnt);

    ATLAS_TRACE_MPI(IRECEIVE) {
        for (idx_t from : neighbours) {
            buffer_requests.push_back(
                comm.iReceive(recv.buffer.data() + recv.displs[from], recv.counts[from], from, buffer_tag));
        }
    }

    ATLAS_TRACE_MPI(WAIT) {
        for (auto request : buffer_requests) {
            comm.wait(request);
        }
    }
#endif
}
}  // namespace

void increase_halo_interior(BuildHaloHelper& helper) {
    helper.update();
    if (helper.node_to_elem.size() == 0) {
        build_lookup_node2elem(helper.mesh, helper.node_to_elem);
    }

    if (helper.uid2node.size() == 0) {
        build_lookup_uid2node(helper.mesh, helper.uid2node);
    }

    // All buffers needed to move elements and nodes
    BuildHaloHelper::Buffers sendmesh(helper.mesh);
    BuildHaloHelper::Buffers recvmesh(helper.mesh);

    // 1) Find boundary nodes of this partition:

    accumulate_partition_bdry_nodes(helper.mesh, helper.halosize, helper.bdry_nodes);
    const std::vector<idx_t>& bdry_nodes = helper.bdry_nodes;
    const idx_t nb_bdry_nodes            = static_cast<idx_t>(bdry_nodes.size());

    // 2) Communicate uid of these boundary nodes to other partitions

    std::vector<uid_t> send_bdry_nodes_uid(bdry_nodes.size());
    for (idx_t jnode = 0; jnode < nb_bdry_nodes; ++jnode) {
        send_bdry_nodes_uid[jnode] = helper.compute_uid(bdry_nodes[jnode]);
    }

    idx_t mpi_size = mpi::size();
    atlas::mpi::Buffer<uid_t, 1> recv_bdry_nodes_uid_from_parts(mpi_size);

    gather_bdry_nodes(helper, send_bdry_nodes_uid, recv_bdry_nodes_uid_from_parts);

    {
        runtime::trace::Barriers set_barriers(false);
        runtime::trace::Logging set_logging(false);
#ifndef ATLAS_103
    /* deprecated */
    atlas_omp_parallel_for (idx_t jpart = 0; jpart < mpi_size; ++jpart)
#else
    const Mesh::PartitionGraph::Neighbours neighbours = helper.mesh.nearestNeighbourPartitions();
    for (idx_t jpart : neighbours)
#endif
    {
        // 3) Find elements and nodes completing these elements in
        //    other tasks that have my nodes through its UID

        mpi::BufferView<uid_t> recv_bdry_nodes_uid = recv_bdry_nodes_uid_from_parts[jpart];

        std::vector<idx_t> found_bdry_elems;
        std::set<uid_t> found_bdry_nodes_uid;

        accumulate_elements(helper.mesh, recv_bdry_nodes_uid, helper.uid2node, helper.node_to_elem, found_bdry_elems,
                            found_bdry_nodes_uid);

        // 4) Fill node and element buffers to send back
        helper.fill_sendbuffer(sendmesh, found_bdry_nodes_uid, found_bdry_elems, jpart);
    }
    }

    // 5) Now communicate all buffers
    helper.all_to_all(sendmesh, recvmesh);

// 6) Adapt mesh
#ifdef DEBUG_OUTPUT
    Log::debug() << "increase_halo_interior -- recv: \n" << recvmesh << std::endl;
#endif
    helper.add_buffers(recvmesh);
}

class PeriodicPoints {
public:
    PeriodicPoints(Mesh& mesh, int flag, idx_t N): flags_(array::make_view<int, 1>(mesh.nodes().flags())) {
        flag_ = flag;
        N_    = N;
    }

    bool operator()(int j) const {
        if (j >= N_) {
            return false;
        }
        if (Topology::check(flags_(j), flag_)) {
            return true;
        }
        return false;
    }

private:
    int N_;
    int flag_;
    array::ArrayView<int, 1> flags_;

    friend std::ostream& operator<<(std::ostream& os, const PeriodicPoints& periodic_points) {
        os << "[";
        for (idx_t j = 0; j < periodic_points.flags_.shape(0); ++j) {
            if (periodic_points(j)) {
                os << " " << j + 1;
            }
        }
        os << " ]";
        return os;
    }
};

void increase_halo_periodic(BuildHaloHelper& helper, const PeriodicPoints& periodic_points,
                            const util::PeriodicTransform& transform, int newflags) {
    helper.update();
    // if (helper.node_to_elem.size() == 0 ) !!! NOT ALLOWED !!! (atlas_test_halo
    // will fail)
    build_lookup_node2elem(helper.mesh, helper.node_to_elem);

    // if( helper.uid2node.size() == 0 ) !!! NOT ALLOWED !!! (atlas_test_halo will
    // fail)
    build_lookup_uid2node(helper.mesh, helper.uid2node);

    // All buffers needed to move elements and nodes
    BuildHaloHelper::Buffers sendmesh(helper.mesh);
    BuildHaloHelper::Buffers recvmesh(helper.mesh);

    // 1) Find boundary nodes of this partition:

    if (!helper.bdry_nodes.size()) {
        accumulate_partition_bdry_nodes(helper.mesh, helper.halosize, helper.bdry_nodes);
    }

    std::vector<idx_t> bdry_nodes = filter_nodes(helper.bdry_nodes, periodic_points);
    const idx_t nb_bdry_nodes     = static_cast<idx_t>(bdry_nodes.size());

    // 2) Compute transformed uid of these boundary nodes and send to other
    // partitions

    std::vector<uid_t> send_bdry_nodes_uid(nb_bdry_nodes);
    for (idx_t jnode = 0; jnode < nb_bdry_nodes; ++jnode) {
        double crd[] = {helper.xy(bdry_nodes[jnode], XX), helper.xy(bdry_nodes[jnode], YY)};
        transform(crd, +1);
        // Log::info() << " crd  " << crd[0] << "  " << crd[1] <<  "       uid " <<
        // util::unique_lonlat(crd) << std::endl;
        send_bdry_nodes_uid[jnode] = util::unique_lonlat(crd);
    }

    idx_t mpi_size = mpi::size();
    atlas::mpi::Buffer<uid_t, 1> recv_bdry_nodes_uid_from_parts(mpi_size);

    gather_bdry_nodes(helper, send_bdry_nodes_uid, recv_bdry_nodes_uid_from_parts,
                      /* periodic = */ true);

    {
        runtime::trace::Barriers set_barriers(false);
        runtime::trace::Logging set_logging(false);
#ifndef ATLAS_103
    /* deprecated */
    atlas_omp_parallel_for (idx_t jpart = 0; jpart < mpi_size; ++jpart)
#else
    Mesh::PartitionGraph::Neighbours neighbours = helper.mesh.nearestNeighbourPartitions();
    // add own rank to neighbours to allow periodicity with self (pole caps)
    idx_t rank = mpi::rank();
    neighbours.insert(std::upper_bound(neighbours.begin(), neighbours.end(), rank), rank);
    for (idx_t jpart : neighbours)
#endif
    {

        // 3) Find elements and nodes completing these elements in
        //    other tasks that have my nodes through its UID

        atlas::mpi::BufferView<uid_t> recv_bdry_nodes_uid = recv_bdry_nodes_uid_from_parts[jpart];

        std::vector<idx_t> found_bdry_elems;
        std::set<uid_t> found_bdry_nodes_uid;

        accumulate_elements(helper.mesh, recv_bdry_nodes_uid, helper.uid2node, helper.node_to_elem, found_bdry_elems,
                            found_bdry_nodes_uid);

        // 4) Fill node and element buffers to send back
        helper.fill_sendbuffer(sendmesh, found_bdry_nodes_uid, found_bdry_elems, transform, newflags, jpart);
    }
    }

    // 5) Now communicate all buffers
    helper.all_to_all(sendmesh, recvmesh);

// 6) Adapt mesh
#ifdef DEBUG_OUTPUT
    Log::debug() << "increase_halo_periodic -- recv: \n" << recvmesh << std::endl;
#endif
    helper.add_buffers(recvmesh);
}

BuildHalo::BuildHalo(Mesh& mesh): mesh_(mesh), periodic_cells_local_index_(mesh.cells().nb_types()) {}

void BuildHalo::operator()(int nb_elems) {
    ATLAS_TRACE("BuildHalo");

    int halo = 0;
    mesh_.metadata().get("halo", halo);

    // Locked halos must be set at mesh generation.
    bool haloLocked = false;
    if (mesh_.metadata().get("halo_locked", haloLocked)) {
        if (haloLocked && halo < nb_elems) {
            const auto errMsg =
                "Error: Halo size locked. Please set MeshGenerator "
                "parameter \"halo\" to at least " +
                std::to_string(nb_elems) +
                ".\n"
                "Current halo size is " +
                std::to_string(halo) + ".\n";
            throw_Exception(errMsg, Here());
        }
    }

    if (halo == nb_elems) {
        return;
    }

    ATLAS_TRACE("Increasing mesh halo");

    for (int jhalo = halo; jhalo < nb_elems; ++jhalo) {
        Log::debug() << "Increase halo " << jhalo + 1 << std::endl;
        idx_t nb_nodes_before_halo_increase = mesh_.nodes().size();

        BuildHaloHelper helper(*this, mesh_);

        ATLAS_TRACE_SCOPE("increase_halo_interior") { increase_halo_interior(helper); }

        PeriodicPoints westpts(mesh_, Topology::PERIODIC | Topology::WEST, nb_nodes_before_halo_increase);

#ifdef DEBUG_OUTPUT
        Log::debug() << "  periodic west : " << westpts << std::endl;
#endif
        ATLAS_TRACE_SCOPE("increase_halo_periodic West") {
            increase_halo_periodic(helper, westpts, WestEast(), Topology::PERIODIC | Topology::WEST | Topology::GHOST);
        }

        PeriodicPoints eastpts(mesh_, Topology::PERIODIC | Topology::EAST, nb_nodes_before_halo_increase);

#ifdef DEBUG_OUTPUT
        Log::debug() << "  periodic east : " << eastpts << std::endl;
#endif
        ATLAS_TRACE_SCOPE("increase_halo_periodic East") {
            increase_halo_periodic(helper, eastpts, EastWest(), Topology::PERIODIC | Topology::EAST | Topology::GHOST);
        }

        for (idx_t p : helper.status.new_periodic_ghost_points) {
            periodic_points_local_index_.push_back(p);
        }
        for (int t = 0; t < mesh_.cells().nb_types(); ++t) {
            for (idx_t p : helper.status.new_periodic_ghost_cells[t]) {
                periodic_cells_local_index_[t].push_back(p);
            }
        }

        {
            std::stringstream ss;
            ss << "nb_nodes_including_halo[" << jhalo + 1 << "]";
            mesh_.metadata().set(ss.str(), mesh_.nodes().size());
        }

        for (int t = 0; t < mesh_.cells().nb_types(); ++t) {
            std::stringstream ss;
            ss << "nb_cells_including_halo[" << t << "][" << jhalo + 1 << "]";
            mesh_.metadata().set(ss.str(), mesh_.cells().elements(t).size());
        }
        mesh_.metadata().set("halo", jhalo + 1);
        mesh_.nodes().global_index().metadata().set("human_readable", false);
        mesh_.cells().global_index().metadata().set("human_readable", false);

#ifdef DEBUG_OUTPUT
        output::Gmsh gmsh2d("build-halo-mesh2d.msh", util::Config("ghost", true)("coordinates", "xy"));
        output::Gmsh gmsh3d("build-halo-mesh3d.msh", util::Config("ghost", true)("coordinates", "xyz"));
        renumber_nodes_glb_idx(mesh_.nodes());
        BuildXYZField("xyz", true)(mesh_.nodes());
        mesh_.metadata().set("halo", jhalo + 1);
        gmsh2d.write(mesh_);
        gmsh3d.write(mesh_);
#endif
    }

    make_nodes_global_index_human_readable(*this, mesh_.nodes(),
                                           /*do_all*/ false);

    make_cells_global_index_human_readable(*this, mesh_.cells(),
                                           /*do_all*/ false);
    //  renumber_nodes_glb_idx (mesh_.nodes());
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

void atlas__build_halo(Mesh::Implementation* mesh, int nb_elems) {
    ATLAS_ASSERT(mesh != nullptr, "Cannot access uninitialised atlas_Mesh");
    Mesh m(mesh);
    build_halo(m, nb_elems);
}

// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
