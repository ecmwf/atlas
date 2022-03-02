/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/actions/BuildEdges.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>

#include "eckit/types/FloatCompare.h"

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/domain.h"
#include "atlas/field/Field.h"
#include "atlas/grid/StructuredGrid.h"
#include "atlas/library/config.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/detail/AccumulateFacets.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/MicroDeg.h"
#include "atlas/util/Unique.h"

using atlas::mesh::detail::accumulate_facets_ordered_by_halo;
using Topology = atlas::mesh::Nodes::Topology;
using atlas::util::microdeg;
using atlas::util::UniqueLonLat;

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

namespace {  // anonymous
struct Sort {
    Sort() = default;
    Sort(gidx_t gid, idx_t idx) {
        g = gid;
        i = idx;
    }
    gidx_t g;
    idx_t i;
    bool operator<(const Sort& other) const { return (g < other.g); }
};
}  // anonymous namespace

void build_element_to_edge_connectivity(Mesh& mesh) {
    ATLAS_TRACE();
    mesh::HybridElements::Connectivity& cell_edge_connectivity = mesh.cells().edge_connectivity();
    cell_edge_connectivity.clear();

    // Allocate cell_edge_connectivity
    for (idx_t t = 0; t < mesh.cells().nb_types(); ++t) {
        idx_t nb_elements       = mesh.cells().elements(t).size();
        idx_t nb_edges_per_elem = mesh.cells().element_type(t).nb_edges();
        std::vector<idx_t> init(mesh.cells().elements(t).size() * nb_edges_per_elem,
                                cell_edge_connectivity.missing_value());
        cell_edge_connectivity.add(nb_elements, nb_edges_per_elem, init.data());
    }

    idx_t nb_edges                                                   = mesh.edges().size();
    const mesh::HybridElements::Connectivity& edge_cell_connectivity = mesh.edges().cell_connectivity();
    const mesh::HybridElements::Connectivity& edge_node_connectivity = mesh.edges().node_connectivity();

    auto edge_flags   = array::make_view<int, 1>(mesh.edges().flags());
    auto is_pole_edge = [&](idx_t e) { return Topology::check(edge_flags(e), Topology::POLE); };

    // Sort edges for bit-reproducibility
    std::vector<Sort> edge_sort;
    edge_sort.reserve(nb_edges);
    {
        UniqueLonLat compute_uid(mesh);

        for (idx_t jedge = 0; jedge < nb_edges; ++jedge) {
            edge_sort.emplace_back(Sort(compute_uid(edge_node_connectivity.row(jedge)), jedge));
        }

        std::sort(edge_sort.data(), edge_sort.data() + nb_edges);
    }

    // Fill in cell_edge_connectivity
    std::vector<idx_t> edge_cnt(mesh.cells().size());
    for (idx_t jedge = 0; jedge < nb_edges; ++jedge) {
        int iedge = edge_sort[jedge].i;
        for (idx_t j = 0; j < 2; ++j) {
            idx_t elem = edge_cell_connectivity(iedge, j);

            if (elem != edge_cell_connectivity.missing_value()) {
                ATLAS_ASSERT(edge_cnt[elem] < cell_edge_connectivity.cols(elem));
                cell_edge_connectivity.set(elem, edge_cnt[elem]++, iedge);
            }
            else {
                if (not is_pole_edge(iedge)) {
                    if (j == 0) {
                        auto node_gidx = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
                        std::stringstream ss;
                        ss << "Edge [" << node_gidx(edge_node_connectivity(jedge, 0)) << ", "
                           << node_gidx(edge_node_connectivity(jedge, 1)) << "] "
                           << "has no element connected.";
                        Log::error() << ss.str() << std::endl;
                        throw_Exception(ss.str(), Here());
                    }
                }
            }
        }
    }


    // Verify that all edges have been found
    auto field_flags = array::make_view<int, 1>(mesh.cells().flags());
    auto patch       = [&field_flags](idx_t e) {
        using Topology = atlas::mesh::Nodes::Topology;
        return Topology::check(field_flags(e), Topology::PATCH);
    };

    for (idx_t jcell = 0; jcell < mesh.cells().size(); ++jcell) {
        if (patch(jcell)) {
            continue;
        }
        for (idx_t jcol = 0; jcol < cell_edge_connectivity.cols(jcell); ++jcol) {
            if (cell_edge_connectivity(jcell, jcol) == cell_edge_connectivity.missing_value()) {
                const array::ArrayView<gidx_t, 1> gidx = array::make_view<gidx_t, 1>(mesh.nodes().global_index());
                std::stringstream msg;
                msg << "Could not find edge " << jcol << " for " << mesh.cells().name(jcell) << " elem " << jcell
                    << " with nodes ( ";
                for (idx_t jnode = 0; jnode < mesh.cells().node_connectivity().cols(jcell); ++jnode) {
                    msg << gidx(mesh.cells().node_connectivity()(jcell, jnode)) << " ";
                }
                msg << ")";
                throw_Exception(msg.str(), Here());
            }
        }
    }
}

void build_node_to_edge_connectivity(Mesh& mesh) {
    mesh::Nodes& nodes   = mesh.nodes();
    const idx_t nb_edges = mesh.edges().size();

    mesh::Nodes::Connectivity& node_to_edge = nodes.edge_connectivity();
    node_to_edge.clear();

    const mesh::HybridElements::Connectivity& edge_node_connectivity = mesh.edges().node_connectivity();

    std::vector<idx_t> to_edge_size(nodes.size(), 0);
    for (idx_t jedge = 0; jedge < nb_edges; ++jedge) {
        for (idx_t j = 0; j < 2; ++j) {
            ++to_edge_size[edge_node_connectivity(jedge, j)];
        }
    }

    node_to_edge.add(nodes.size(), to_edge_size.data());

    for (idx_t jnode = 0; jnode < nodes.size(); ++jnode) {
        to_edge_size[jnode] = 0;
    }

    UniqueLonLat compute_uid(mesh);
    std::vector<Sort> edge_sort(nb_edges);
    for (idx_t jedge = 0; jedge < nb_edges; ++jedge) {
        edge_sort[jedge] = Sort(compute_uid(edge_node_connectivity.row(jedge)), jedge);
    }
    std::stable_sort(edge_sort.data(), edge_sort.data() + nb_edges);

    for (idx_t jedge = 0; jedge < nb_edges; ++jedge) {
        idx_t iedge = edge_sort[jedge].i;
        ATLAS_ASSERT(iedge < nb_edges);
        for (idx_t j = 0; j < 2; ++j) {
            idx_t node = edge_node_connectivity(iedge, j);
            node_to_edge.set(node, to_edge_size[node]++, iedge);
        }
    }
}

class AccumulatePoleEdges {
    enum
    {
        NORTH = 0,
        SOUTH = 1
    };
    const array::ArrayView<double, 2> xy;
    const array::ArrayView<int, 1> flags;
    const array::ArrayView<int, 1> part;
    const array::ArrayView<int, 1> halo;
    const idx_t nb_nodes;
    std::vector<std::set<int>> pole_nodes;

public:
    AccumulatePoleEdges(mesh::Nodes& nodes):
        xy(array::make_view<double, 2>(nodes.xy())),
        flags(array::make_view<int, 1>(nodes.flags())),
        part(array::make_view<int, 1>(nodes.partition())),
        halo(array::make_view<int, 1>(nodes.halo())),
        nb_nodes(nodes.size()),
        pole_nodes(2) {
        double min[2], max[2];
        min[XX] = std::numeric_limits<double>::max();
        min[YY] = std::numeric_limits<double>::max();
        max[XX] = -std::numeric_limits<double>::max();
        max[YY] = -std::numeric_limits<double>::max();
        for (idx_t node = 0; node < nb_nodes; ++node) {
            min[XX] = std::min(min[XX], xy(node, XX));
            min[YY] = std::min(min[YY], xy(node, YY));
            max[XX] = std::max(max[XX], xy(node, XX));
            max[YY] = std::max(max[YY], xy(node, YY));
        }

        ATLAS_TRACE_MPI(ALLREDUCE) {
            mpi::comm().allReduceInPlace(min, 2, eckit::mpi::min());
            mpi::comm().allReduceInPlace(max, 2, eckit::mpi::max());
        }

        double tol = 1e-6;

        // Collect all nodes closest to poles
        for (idx_t node = 0; node < nb_nodes; ++node) {
            if (std::abs(xy(node, YY) - max[YY]) < tol) {
                pole_nodes[NORTH].insert(node);
            }
            else if (std::abs(xy(node, YY) - min[YY]) < tol) {
                pole_nodes[SOUTH].insert(node);
            }
        }

        // Sanity check
        {
            for (idx_t NS = 0; NS < 2; ++NS) {
                int npart = -1;
                for (std::set<int>::iterator it = pole_nodes[NS].begin(); it != pole_nodes[NS].end(); ++it) {
                    int node = *it;
                    if (npart == -1) {
                        npart = part(node);
                    }
                    else if (part(node) != npart) {
                        // Not implemented yet, when pole-lattitude is split.
                        std::stringstream msg;
                        msg << "Split pole-latitude is not supported yet...  node " << node << "[p" << part(node)
                            << "] should belong to part " << npart;
                        throw_NotImplemented(msg.str(), Here());
                    }
                }
            }
        }
    }
    void compute_pole_edges(int _halo, std::vector<idx_t>& pole_edge_nodes, idx_t& nb_pole_edges) {
        // Create connections over the poles and store in pole_edge_nodes
        nb_pole_edges = 0;
        for (idx_t NS = 0; NS < 2; ++NS) {
            for (std::set<int>::iterator it = pole_nodes[NS].begin(); it != pole_nodes[NS].end(); ++it) {
                int node = *it;
                if (!Topology::check(flags(node), Topology::PERIODIC | Topology::GHOST)) {
                    int x2 = microdeg(xy(node, XX) + 180.);
                    for (std::set<int>::iterator itr = pole_nodes[NS].begin(); itr != pole_nodes[NS].end(); ++itr) {
                        int other_node = *itr;
                        if (microdeg(xy(other_node, XX)) == x2) {
                            if (!Topology::check(flags(other_node), Topology::PERIODIC)) {
                                if (halo(node) == _halo && halo(other_node) == _halo) {
                                    pole_edge_nodes.push_back(node);
                                    pole_edge_nodes.push_back(other_node);
                                    ++nb_pole_edges;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};


struct ComputeUniquePoleEdgeIndex {
    // Already assumes that the edges cross the pole

    ComputeUniquePoleEdgeIndex(const mesh::Nodes& nodes): xy(array::make_view<double, 2>(nodes.xy())) {}

    gidx_t operator()(const mesh::Connectivity::Row& edge_nodes) const {
        double centroid[2];
        centroid[XX] = 0.;
        centroid[YY] = 0.;
        for (idx_t jnode = 0; jnode < 2; ++jnode) {
            centroid[XX] += xy(edge_nodes(jnode), XX);
            centroid[YY] += xy(edge_nodes(jnode), YY);
        }
        centroid[XX] /= 2.;
        centroid[YY] /= 2.;
        if (centroid[YY] > 0) {
            centroid[YY] = 90.;
        }
        else {
            centroid[YY] = -90.;
        }
        /// FIXME make this into `util::unique_lonlat(centroid)` but this causes
        /// weird parallel behavior
        return util::detail::unique32(microdeg(centroid[XX]), microdeg(centroid[YY]));
    }

    array::ArrayView<const double, 2> xy;
};

void build_edges(Mesh& mesh) {
    build_edges(mesh, util::NoConfig());
}

void build_edges(Mesh& mesh, const eckit::Configuration& config) {
    ATLAS_TRACE("BuildEdges");

    int mesh_halo(0);
    mesh.metadata().get("halo", mesh_halo);

    if (mesh.metadata().has("built_edges_for_halo")) {
        int edges_halo = mesh.metadata().getInt("built_edges_for_halo");
        if (edges_halo == mesh_halo) {
            // Nothing to be done here
            return;
        }
    }

    bool pole_edges{false};
    if (StructuredGrid grid = mesh.grid()) {
        if (RectangularDomain domain = grid.domain()) {
            if (domain.global()) {
                double ymax = std::max(std::abs(grid.y().front()), std::abs(grid.y().back()));
                pole_edges  = not eckit::types::is_approximately_equal(ymax, domain.ymax());
            }
        }
    }
    config.get("pole_edges", pole_edges);

    bool sort_edges{false};
    config.get("sort_edges", sort_edges);


    mesh::Nodes& nodes = mesh.nodes();
    auto node_part     = array::make_view<int, 1>(nodes.partition());

    idx_t nb_nodes = nodes.size();

    mesh.edges().clear();

    idx_t edge_start{0};
    idx_t edge_end{0};

    // storage for edge-to-node-connectivity shape=(nb_edges,2)
    std::vector<idx_t> edge_nodes_data;
    std::vector<idx_t> edge_to_elem_data;
    std::vector<idx_t> edge_halo_offsets;
    idx_t nb_edges;
    idx_t nb_inner_edges;
    idx_t missing_value;

    accumulate_facets_ordered_by_halo(mesh.cells(), mesh.nodes(), edge_nodes_data, edge_to_elem_data, nb_edges,
                                      nb_inner_edges, missing_value, edge_halo_offsets);

    std::vector<idx_t> sorted_edge_nodes_data;
    std::vector<idx_t> sorted_edge_to_elem_data;

    for (int halo = 0; halo <= mesh_halo; ++halo) {
        edge_start = edge_end;
        edge_end += (edge_halo_offsets[halo + 1] - edge_halo_offsets[halo]);

        if (/*sort edges based on lowest node local index = */ sort_edges) {
            if (sorted_edge_nodes_data.empty()) {
                sorted_edge_nodes_data.resize(edge_nodes_data.size());
            }
            if (sorted_edge_to_elem_data.empty()) {
                sorted_edge_to_elem_data.resize(edge_to_elem_data.size());
            }
            std::vector<std::pair<idx_t, idx_t>> sorted_edges_by_lowest_node_index;
            sorted_edges_by_lowest_node_index.reserve(edge_end - edge_start);
            for (idx_t e = edge_start; e < edge_end; ++e) {
                const idx_t iedge     = edge_halo_offsets[halo] + (e - edge_start);
                idx_t lowest_node_idx = std::min(edge_nodes_data.at(2 * iedge + 0), edge_nodes_data.at(2 * iedge + 1));
                sorted_edges_by_lowest_node_index.emplace_back(lowest_node_idx, e);
            }
            std::sort(sorted_edges_by_lowest_node_index.begin(), sorted_edges_by_lowest_node_index.end());
            for (idx_t e = edge_start; e < edge_end; ++e) {
                const idx_t iedge = edge_halo_offsets[halo] + (e - edge_start);
                const idx_t sedge =
                    edge_halo_offsets[halo] + (sorted_edges_by_lowest_node_index[e - edge_start].second - edge_start);
                sorted_edge_nodes_data[2 * iedge + 0]   = edge_nodes_data[2 * sedge + 0];
                sorted_edge_nodes_data[2 * iedge + 1]   = edge_nodes_data[2 * sedge + 1];
                sorted_edge_to_elem_data[2 * iedge + 0] = edge_to_elem_data[2 * sedge + 0];
                sorted_edge_to_elem_data[2 * iedge + 1] = edge_to_elem_data[2 * sedge + 1];
            }

            for (idx_t e = edge_start; e < edge_end; ++e) {
                const idx_t iedge                = edge_halo_offsets[halo] + (e - edge_start);
                edge_nodes_data[2 * iedge + 0]   = sorted_edge_nodes_data[2 * iedge + 0];
                edge_nodes_data[2 * iedge + 1]   = sorted_edge_nodes_data[2 * iedge + 1];
                edge_to_elem_data[2 * iedge + 0] = sorted_edge_to_elem_data[2 * iedge + 0];
                edge_to_elem_data[2 * iedge + 1] = sorted_edge_to_elem_data[2 * iedge + 1];
            }
        }

        // Build edges
        mesh.edges().add(new mesh::temporary::Line(), (edge_end - edge_start),
                         edge_nodes_data.data() + edge_halo_offsets[halo] * 2);
        auto& edge_nodes       = mesh.edges().node_connectivity();
        const auto& cell_nodes = mesh.cells().node_connectivity();

        UniqueLonLat compute_uid(mesh);

        auto edge_ridx    = array::make_indexview<idx_t, 1>(mesh.edges().remote_index());
        auto edge_part    = array::make_view<int, 1>(mesh.edges().partition());
        auto edge_glb_idx = array::make_view<gidx_t, 1>(mesh.edges().global_index());
        auto edge_halo    = array::make_view<int, 1>(mesh.edges().halo());
        auto edge_flags   = array::make_view<int, 1>(mesh.edges().flags());

        ATLAS_ASSERT(cell_nodes.missing_value() == missing_value);
        for (idx_t edge = edge_start; edge < edge_end; ++edge) {
            const idx_t iedge = edge_halo_offsets[halo] + (edge - edge_start);
            const int ip1     = edge_nodes(edge, 0);
            const int ip2     = edge_nodes(edge, 1);
            if (compute_uid(ip1) > compute_uid(ip2)) {
                idx_t swapped[2] = {ip2, ip1};
                edge_nodes.set(edge, swapped);
            }

            ATLAS_ASSERT(idx_t(edge_nodes(edge, 0)) < nb_nodes);
            ATLAS_ASSERT(idx_t(edge_nodes(edge, 1)) < nb_nodes);
            edge_glb_idx(edge) = compute_uid(edge_nodes.row(edge));
            edge_part(edge)    = std::min(node_part(edge_nodes(edge, 0)), node_part(edge_nodes(edge, 1)));
            edge_ridx(edge)    = edge;
            edge_halo(edge)    = halo;
            edge_flags(edge)   = 0;

            const idx_t e1 = edge_to_elem_data[2 * iedge + 0];
            const idx_t e2 = edge_to_elem_data[2 * iedge + 1];

            ATLAS_ASSERT(e1 != cell_nodes.missing_value());
            if (e2 == cell_nodes.missing_value()) {
                // do nothing
            }
            else if (compute_uid(cell_nodes.row(e1)) > compute_uid(cell_nodes.row(e2))) {
                edge_to_elem_data[iedge * 2 + 0] = e2;
                edge_to_elem_data[iedge * 2 + 1] = e1;
            }
        }

        mesh.edges().cell_connectivity().add((edge_end - edge_start), 2,
                                             edge_to_elem_data.data() + edge_halo_offsets[halo] * 2);

        if (pole_edges) {
            auto pole_edge_accumulator = std::make_shared<AccumulatePoleEdges>(nodes);

            idx_t nb_pole_edges;
            std::vector<idx_t> pole_edge_nodes;

            pole_edge_accumulator->compute_pole_edges(halo, pole_edge_nodes, nb_pole_edges);

            if (nb_pole_edges) {
                edge_start = edge_end;
                edge_end += nb_pole_edges;

                mesh.edges().add(new mesh::temporary::Line(), nb_pole_edges, pole_edge_nodes.data());

                auto edge_ridx    = array::make_indexview<idx_t, 1>(mesh.edges().remote_index());
                auto edge_part    = array::make_view<int, 1>(mesh.edges().partition());
                auto edge_glb_idx = array::make_view<gidx_t, 1>(mesh.edges().global_index());
                auto edge_halo    = array::make_view<int, 1>(mesh.edges().halo());
                auto edge_flags   = array::make_view<int, 1>(mesh.edges().flags());

                auto set_pole_edge = [&edge_flags](idx_t e) { Topology::set(edge_flags(e), Topology::POLE); };

                auto& edge_nodes = mesh.edges().node_connectivity();

                mesh.edges().cell_connectivity().add(nb_pole_edges, 2);

                idx_t cnt = 0;
                ComputeUniquePoleEdgeIndex compute_uid(nodes);
                for (idx_t edge = edge_start; edge < edge_end; ++edge) {
                    idx_t ip1 = pole_edge_nodes[cnt++];
                    idx_t ip2 = pole_edge_nodes[cnt++];
                    std::array<idx_t, 2> enodes{ip1, ip2};
                    edge_nodes.set(edge, enodes.data());
                    edge_glb_idx(edge) = compute_uid(edge_nodes.row(edge));
                    edge_part(edge)    = std::min(node_part(edge_nodes(edge, 0)), node_part(edge_nodes(edge, 1)));
                    edge_ridx(edge)    = edge;
                    edge_halo(edge)    = halo;
                    set_pole_edge(edge);
                }
            }
        }
    }

    mesh.edges().metadata().set("pole_edges", pole_edges);


    build_element_to_edge_connectivity(mesh);

    mesh::HybridElements::Connectivity& cell_edges = mesh.cells().edge_connectivity();
    auto cell_halo                                 = array::make_view<int, 1>(mesh.cells().halo());
    auto cell_flags                                = array::make_view<int, 1>(mesh.cells().flags());
    auto cell_patch                                = [&cell_flags](idx_t e) {
        using Topology = atlas::mesh::Nodes::Topology;
        return Topology::check(cell_flags(e), Topology::PATCH);
    };
    auto edge_halo = array::make_view<int, 1>(mesh.edges().halo());
    int max_halo   = 0;
    for (idx_t jcell = 0; jcell < mesh.cells().size(); ++jcell) {
        if (not cell_patch(jcell)) {
            int halo = cell_halo(jcell);
            max_halo = std::max(halo, max_halo);
            for (idx_t jedge = 0; jedge < cell_edges.cols(jcell); ++jedge) {
                auto iedge = cell_edges(jcell, jedge);
                ATLAS_ASSERT(edge_halo(iedge) <= halo);
            }
        }
    }

    std::vector<int> nb_edges_including_halo(max_halo + 1);

    {
        int nb_edges = mesh.edges().size();
        for (int jedge = 0; jedge < nb_edges; ++jedge) {
            nb_edges_including_halo[edge_halo(jedge)] = jedge + 1;
            if (jedge > 0) {
                ATLAS_ASSERT(edge_halo(jedge) >= edge_halo(jedge - 1));
            }
        }
    }

    for (int i = 0; i <= max_halo; ++i) {
        if (i > 0) {
            ATLAS_ASSERT(nb_edges_including_halo[i] > nb_edges_including_halo[i - 1]);
        }
        std::stringstream ss;
        ss << "nb_edges_including_halo[" << i << "]";
        mesh.metadata().set(ss.str(), nb_edges_including_halo[i]);
    }

    mesh.metadata().set("built_edges_for_halo", mesh_halo);

    // Backwards compatibility for code that reads "is_pole_edge" field instead of checking the flags, only Fortran would do it
    {
        if (pole_edges) {
            if (!mesh.edges().has_field("is_pole_edge")) {
                mesh.edges().add(
                    Field("is_pole_edge", array::make_datatype<int>(), array::make_shape(mesh.edges().size())));
            }
            auto edge_flags   = array::make_view<int, 1>(mesh.edges().flags());
            auto is_pole_edge = array::make_view<int, 1>(mesh.edges().field("is_pole_edge"));
            int nb_edges      = mesh.edges().size();
            for (int jedge = 0; jedge < nb_edges; ++jedge) {
                is_pole_edge(jedge) = Topology::check(edge_flags(jedge), Topology::POLE);
            }
        }
    }
}

void build_pole_edges(Mesh&) {
    ATLAS_TRACE();
    Log::info() << "ATLAS_WARNING: Deprecation warning: build_pole_edges is no longer required.\n"
                << "It is automatically inferred within atlas_build_edges" << std::endl;
    Log::info() << "The 'build_pole_edges' function will be removed in a future version" << std::endl;
}

//----------------------------------------------------------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {
void atlas__build_edges(Mesh::Implementation* mesh) {
    ATLAS_ASSERT(mesh != nullptr, "Cannot access uninitialised atlas_Mesh");
    Mesh m(mesh);
    build_edges(m);
}
void atlas__build_pole_edges(Mesh::Implementation*) {
    Log::info() << "ATLAS_WARNING: Deprecation warning: atlas_build_pole_edges is no longer required.\n"
                << "It is automatically inferred within atlas_build_edges" << std::endl;
    Log::info() << "The 'atlas_build_pole_edges' function will be removed in a future version" << std::endl;
}
void atlas__build_node_to_edge_connectivity(Mesh::Implementation* mesh) {
    ATLAS_ASSERT(mesh != nullptr, "Cannot access uninitialised atlas_Mesh");
    Mesh m(mesh);
    build_node_to_edge_connectivity(m);
}
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
