/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>
#include <iomanip>
#include <numeric>

#include "eckit/filesystem/PathName.h"

#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/IsGhostNode.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/WriteLoadBalanceReport.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/library/FloatingPointExceptions.h"

using atlas::mesh::IsGhostNode;

namespace atlas {
namespace mesh {
namespace actions {

void write_load_balance_report(const Mesh& mesh, const std::string& filename) {
    std::ofstream ofs;
    if (mpi::rank() == 0) {
        eckit::PathName path(filename);
        ofs.open(path.localPath(), std::ofstream::out);
    }

    write_load_balance_report(mesh, ofs);

    if (mpi::rank() == 0) {
        ofs.close();
    }
}

void write_load_balance_report(const Mesh& mesh, std::ostream& ofs) {
    idx_t npart = mpi::size();
    idx_t root  = 0;

    std::vector<idx_t> nb_total_nodes(npart, 0);
    std::vector<idx_t> nb_owned_nodes(npart, 0);
    std::vector<idx_t> nb_ghost_nodes(npart, 0);
    std::vector<double> ghost_ratio_nodes(npart, 0);

    std::vector<idx_t> nb_total_edges(npart, 0);
    std::vector<idx_t> nb_owned_edges(npart, 0);
    std::vector<idx_t> nb_ghost_edges(npart, 0);
    std::vector<double> nb_ghost_ratio_edges(npart, 0);

    {
        const mesh::Nodes& nodes = mesh.nodes();
        IsGhostNode is_ghost(nodes);
        idx_t nb_nodes = nodes.size();
        idx_t nowned(0);
        idx_t nghost(0);
        for (idx_t n = 0; n < nb_nodes; ++n) {
            if (is_ghost(n)) {
                ++nghost;
            }
            else {
                ++nowned;
            }
        }

        /// @note this could be improved by packing the 3 integers in a vector, and
        /// doing only comm() call

        ATLAS_TRACE_MPI(GATHER) {
            mpi::comm().gather(nb_nodes, nb_total_nodes, root);
            mpi::comm().gather(nowned, nb_owned_nodes, root);
            mpi::comm().gather(nghost, nb_ghost_nodes, root);
        }

        bool disabled_fpe = library::disable_floating_point_exception(FE_INVALID);

        for (idx_t p = 0; p < npart; ++p) {
            if (nb_owned_nodes[p]) {
                ghost_ratio_nodes[p] = static_cast<double>(nb_ghost_nodes[p]) / static_cast<double>(nb_owned_nodes[p]);
            }
            else {
                ghost_ratio_nodes[p] = -1;
            }
        }
        if (disabled_fpe) {
            library::enable_floating_point_exception(FE_INVALID);
        }
    }

    bool has_edges = mesh.edges().size();

    if (has_edges) {
        const mesh::Nodes& nodes = mesh.nodes();
        IsGhostNode is_ghost(nodes);
        const mesh::HybridElements::Connectivity& edge_nodes = mesh.edges().node_connectivity();
        idx_t nb_edges                                       = mesh.edges().size();
        idx_t nowned(0);
        idx_t nghost(0);
        for (idx_t j = 0; j < nb_edges; ++j) {
            if (is_ghost(edge_nodes(j, 0))) {
                ++nghost;
            }
            else {
                ++nowned;
            }
        }

        /// @note this could be improved by packing the 3 integers in a vector, and
        /// doing only comm() call

        mpi::comm().gather(nb_edges, nb_total_edges, root);
        mpi::comm().gather(nowned, nb_owned_edges, root);
        mpi::comm().gather(nghost, nb_ghost_nodes, root);
    }

    if (mpi::rank() == 0) {
        int idt = 10;
        ofs << "# STATISTICS\n";
        ofs << std::setw(1) << "#" << std::setw(5) << "";
        ofs << std::setw(idt) << "nodes";
        ofs << std::setw(idt) << "owned";
        ofs << std::setw(idt) << "ghost";
        ofs << std::setw(idt) << "ratio(%)";
        if (has_edges) {
            ofs << std::setw(idt) << "edges";
            ofs << std::setw(idt) << "oedges";
            ofs << std::setw(idt) << "gedges";
        }
        ofs << "\n";
        ofs << std::setw(6) << "# tot ";
        ofs << std::setw(idt) << std::accumulate(nb_total_nodes.data(), nb_total_nodes.data() + npart, 0);
        ofs << std::setw(idt) << std::accumulate(nb_owned_nodes.data(), nb_owned_nodes.data() + npart, 0);
        ofs << std::setw(idt) << std::accumulate(nb_ghost_nodes.data(), nb_ghost_nodes.data() + npart, 0);
        ofs << std::setw(idt) << "/";
        if (has_edges) {
            ofs << std::setw(idt) << std::accumulate(nb_total_edges.data(), nb_total_edges.data() + npart, 0);
            ofs << std::setw(idt) << std::accumulate(nb_owned_edges.data(), nb_owned_edges.data() + npart, 0);
            ofs << std::setw(idt) << std::accumulate(nb_ghost_edges.data(), nb_ghost_edges.data() + npart, 0);
        }
        ofs << "\n";
        ofs << std::setw(6) << "# max ";
        ofs << std::setw(idt) << *std::max_element(nb_total_nodes.data(), nb_total_nodes.data() + npart);
        ofs << std::setw(idt) << *std::max_element(nb_owned_nodes.data(), nb_owned_nodes.data() + npart);
        ofs << std::setw(idt) << *std::max_element(nb_ghost_nodes.data(), nb_ghost_nodes.data() + npart);
        ofs << std::setw(idt) << std::setw(idt) << std::fixed << std::setprecision(2)
            << *std::max_element(ghost_ratio_nodes.data(), ghost_ratio_nodes.data() + npart) * 100.;
        if (has_edges) {
            ofs << std::setw(idt) << *std::max_element(nb_total_edges.data(), nb_total_edges.data() + npart);
            ofs << std::setw(idt) << *std::max_element(nb_owned_edges.data(), nb_owned_edges.data() + npart);
            ofs << std::setw(idt) << *std::max_element(nb_ghost_edges.data(), nb_ghost_edges.data() + npart);
        }
        ofs << "\n";
        ofs << std::setw(6) << "# min ";
        ofs << std::setw(idt) << *std::min_element(nb_total_nodes.data(), nb_total_nodes.data() + npart);
        ofs << std::setw(idt) << *std::min_element(nb_owned_nodes.data(), nb_owned_nodes.data() + npart);
        ofs << std::setw(idt) << *std::min_element(nb_ghost_nodes.data(), nb_ghost_nodes.data() + npart);
        ofs << std::setw(idt) << std::fixed << std::setprecision(2)
            << *std::min_element(ghost_ratio_nodes.data(), ghost_ratio_nodes.data() + npart) * 100.;
        if (has_edges) {
            ofs << std::setw(idt) << *std::min_element(nb_total_edges.data(), nb_total_edges.data() + npart);
            ofs << std::setw(idt) << *std::min_element(nb_owned_edges.data(), nb_owned_edges.data() + npart);
            ofs << std::setw(idt) << *std::min_element(nb_ghost_edges.data(), nb_ghost_edges.data() + npart);
        }
        ofs << "\n";
        ofs << std::setw(6) << "# avg ";
        ofs << std::setw(idt) << std::accumulate(nb_total_nodes.data(), nb_total_nodes.data() + npart, 0) / npart;
        ofs << std::setw(idt) << std::accumulate(nb_owned_nodes.data(), nb_owned_nodes.data() + npart, 0) / npart;
        ofs << std::setw(idt) << std::accumulate(nb_ghost_nodes.data(), nb_ghost_nodes.data() + npart, 0) / npart;
        ofs << std::setw(idt) << std::fixed << std::setprecision(2)
            << std::accumulate(ghost_ratio_nodes.data(), ghost_ratio_nodes.data() + npart, 0.) /
                   static_cast<double>(npart) * 100.;
        if (has_edges) {
            ofs << std::setw(idt) << std::accumulate(nb_total_edges.data(), nb_total_edges.data() + npart, 0) / npart;
            ofs << std::setw(idt) << std::accumulate(nb_owned_edges.data(), nb_owned_edges.data() + npart, 0) / npart;
            ofs << std::setw(idt) << std::accumulate(nb_ghost_edges.data(), nb_ghost_edges.data() + npart, 0) / npart;
        }
        ofs << "\n";
        ofs << "#----------------------------------------------------\n";
        ofs << "# PER TASK\n";
        ofs << std::setw(6) << "# part";
        ofs << std::setw(idt) << "nodes";
        ofs << std::setw(idt) << "owned";
        ofs << std::setw(idt) << "ghost";
        ofs << std::setw(idt) << "ratio(%)";
        if (has_edges) {
            ofs << std::setw(idt) << "edges";
            ofs << std::setw(idt) << "oedges";
            ofs << std::setw(idt) << "gedges";
        }
        ofs << "\n";
        for (idx_t jpart = 0; jpart < npart; ++jpart) {
            ofs << std::setw(6) << jpart;
            ofs << std::setw(idt) << nb_total_nodes[jpart];
            ofs << std::setw(idt) << nb_owned_nodes[jpart];
            ofs << std::setw(idt) << nb_ghost_nodes[jpart];
            ofs << std::setw(idt) << std::fixed << std::setprecision(2) << ghost_ratio_nodes[jpart] * 100.;
            if (has_edges) {
                ofs << std::setw(idt) << nb_total_edges[jpart];
                ofs << std::setw(idt) << nb_owned_edges[jpart];
                ofs << std::setw(idt) << nb_ghost_edges[jpart];
            }
            ofs << "\n";
        }
    }
}

// ------------------------------------------------------------------

// C wrapper interfaces to C++ routines
void atlas__write_load_balance_report(Mesh::Implementation* mesh, char* filename) {
    ATLAS_ASSERT(mesh != nullptr, "Cannot access uninitialised atlas_Mesh");
    Mesh m(mesh);
    write_load_balance_report(m, std::string(filename));
}

// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
