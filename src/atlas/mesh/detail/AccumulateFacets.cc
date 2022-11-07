/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/mesh/detail/AccumulateFacets.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"

namespace atlas {
namespace mesh {
namespace detail {

void accumulate_facets(const mesh::HybridElements& cells, const mesh::Nodes& nodes,
                       std::vector<idx_t>& facet_nodes_data,  // shape(nb_facets,nb_nodes_per_facet)
                       std::vector<idx_t>& connectivity_facet_to_elem, idx_t& nb_facets, idx_t& nb_inner_facets,
                       idx_t& missing_value) {
    ATLAS_TRACE();
    missing_value = -1;
    std::vector<std::vector<idx_t>> node_to_facet(nodes.size());
    for (auto& facet : node_to_facet) {
        facet.reserve(6);
    }
    nb_facets       = 0;
    nb_inner_facets = 0;
    if (connectivity_facet_to_elem.size() == 0) {
        connectivity_facet_to_elem.reserve(6 * cells.size());
    }
    if (facet_nodes_data.size() == 0) {
        facet_nodes_data.reserve(6 * cells.size());
    }
    for (idx_t t = 0; t < cells.nb_types(); ++t) {
        const mesh::Elements& elements            = cells.elements(t);
        const mesh::BlockConnectivity& elem_nodes = elements.node_connectivity();
        auto elem_flags                           = elements.view<int, 1>(elements.flags());

        auto patch = [&elem_flags](idx_t e) {
            using Topology = atlas::mesh::Nodes::Topology;
            return Topology::check(elem_flags(e), Topology::PATCH);
        };

        idx_t nb_elems          = elements.size();
        idx_t nb_nodes_in_facet = 2;

        std::vector<std::vector<int>> facet_node_numbering;
        idx_t nb_facets_in_elem;
        if (elements.name() == "Pentagon") {
            nb_facets_in_elem = 5;
            facet_node_numbering.resize(nb_facets_in_elem, std::vector<int>(nb_nodes_in_facet));
            facet_node_numbering[0][0] = 0;
            facet_node_numbering[0][1] = 1;
            facet_node_numbering[1][0] = 1;
            facet_node_numbering[1][1] = 2;
            facet_node_numbering[2][0] = 2;
            facet_node_numbering[2][1] = 3;
            facet_node_numbering[3][0] = 3;
            facet_node_numbering[3][1] = 4;
            facet_node_numbering[4][0] = 4;
            facet_node_numbering[4][1] = 0;
        }
        else if (elements.name() == "Quadrilateral") {
            nb_facets_in_elem = 4;
            facet_node_numbering.resize(nb_facets_in_elem, std::vector<int>(nb_nodes_in_facet));
            facet_node_numbering[0][0] = 0;
            facet_node_numbering[0][1] = 1;
            facet_node_numbering[1][0] = 1;
            facet_node_numbering[1][1] = 2;
            facet_node_numbering[2][0] = 2;
            facet_node_numbering[2][1] = 3;
            facet_node_numbering[3][0] = 3;
            facet_node_numbering[3][1] = 0;
        }
        else if (elements.name() == "Triangle") {
            nb_facets_in_elem = 3;
            facet_node_numbering.resize(nb_facets_in_elem, std::vector<int>(nb_nodes_in_facet));
            facet_node_numbering[0][0] = 0;
            facet_node_numbering[0][1] = 1;
            facet_node_numbering[1][0] = 1;
            facet_node_numbering[1][1] = 2;
            facet_node_numbering[2][0] = 2;
            facet_node_numbering[2][1] = 0;
        }
        else {
            throw_Exception(elements.name() + " is not \"Pentagon\", \"Quadrilateral\", or \"Triangle\"", Here());
        }

        std::vector<idx_t> facet_nodes(nb_nodes_in_facet);

        for (idx_t e = 0; e < nb_elems; ++e) {
            if (patch(e)) {
                continue;
            }
            for (idx_t f = 0; f < nb_facets_in_elem; ++f) {
                bool found_face = false;

                for (idx_t jnode = 0; jnode < nb_nodes_in_facet; ++jnode) {
                    facet_nodes[jnode] = elem_nodes(e, facet_node_numbering[f][jnode]);
                }

                idx_t node = facet_nodes[0];
                for (const idx_t face : node_to_facet[node]) {
                    idx_t nb_matched_nodes = 0;
                    if (nb_nodes_in_facet > 1)  // 2D or 3D
                    {
                        for (idx_t jnode = 0; jnode < nb_nodes_in_facet; ++jnode) {
                            idx_t other_node = facet_nodes[jnode];
                            for (const idx_t other_face : node_to_facet[other_node]) {
                                if (other_face == face) {
                                    ++nb_matched_nodes;
                                    break;
                                }
                            }
                        }
                        if (nb_matched_nodes == nb_nodes_in_facet) {
                            connectivity_facet_to_elem[2 * face + 1] = e + elements.begin();
                            ++nb_inner_facets;
                            found_face = true;
                            break;
                        }
                    }
                }

                if (found_face == false) {
                    connectivity_facet_to_elem.emplace_back(elements.begin() + e);
                    // if 2nd element stays missing_value, it is a bdry face
                    connectivity_facet_to_elem.emplace_back(missing_value);
                    for (idx_t n = 0; n < nb_nodes_in_facet; ++n) {
                        node_to_facet[facet_nodes[n]].emplace_back(nb_facets);
                        facet_nodes_data.emplace_back(facet_nodes[n]);
                    }
                    ++nb_facets;
                }
            }
        }
    }
}

void accumulate_facets_in_range(std::vector<array::Range>& range, const mesh::HybridElements& cells,
                                const mesh::Nodes& /*nodes*/,
                                std::vector<idx_t>& facet_nodes_data,  // shape(nb_facets,nb_nodes_per_facet)
                                std::vector<idx_t>& connectivity_facet_to_elem, idx_t& nb_facets,
                                idx_t& nb_inner_facets, idx_t& missing_value,
                                std::vector<std::vector<idx_t>>& node_to_facet) {
    ATLAS_TRACE();
    if (connectivity_facet_to_elem.size() == 0) {
        connectivity_facet_to_elem.reserve(6 * cells.size());
    }
    if (facet_nodes_data.size() == 0) {
        facet_nodes_data.reserve(6 * cells.size());
    }
    for (idx_t t = 0; t < cells.nb_types(); ++t) {
        const mesh::Elements& elements            = cells.elements(t);
        const mesh::BlockConnectivity& elem_nodes = elements.node_connectivity();
        auto elem_flags                           = elements.view<int, 1>(elements.flags());

        auto patch = [&elem_flags](idx_t e) {
            using Topology = atlas::mesh::Nodes::Topology;
            return Topology::check(elem_flags(e), Topology::PATCH);
        };

        idx_t nb_nodes_in_facet = 2;

        std::vector<std::vector<int>> facet_node_numbering;
        idx_t nb_facets_in_elem;
        if (elements.name() == "Pentagon") {
            nb_facets_in_elem = 5;
            facet_node_numbering.resize(nb_facets_in_elem, std::vector<int>(nb_nodes_in_facet));
            facet_node_numbering[0][0] = 0;
            facet_node_numbering[0][1] = 1;
            facet_node_numbering[1][0] = 1;
            facet_node_numbering[1][1] = 2;
            facet_node_numbering[2][0] = 2;
            facet_node_numbering[2][1] = 3;
            facet_node_numbering[3][0] = 3;
            facet_node_numbering[3][1] = 4;
            facet_node_numbering[4][0] = 4;
            facet_node_numbering[4][1] = 0;
        }
        else if (elements.name() == "Quadrilateral") {
            nb_facets_in_elem = 4;
            facet_node_numbering.resize(nb_facets_in_elem, std::vector<int>(nb_nodes_in_facet));
            facet_node_numbering[0][0] = 0;
            facet_node_numbering[0][1] = 1;
            facet_node_numbering[1][0] = 1;
            facet_node_numbering[1][1] = 2;
            facet_node_numbering[2][0] = 2;
            facet_node_numbering[2][1] = 3;
            facet_node_numbering[3][0] = 3;
            facet_node_numbering[3][1] = 0;
        }
        else if (elements.name() == "Triangle") {
            nb_facets_in_elem = 3;
            facet_node_numbering.resize(nb_facets_in_elem, std::vector<int>(nb_nodes_in_facet));
            facet_node_numbering[0][0] = 0;
            facet_node_numbering[0][1] = 1;
            facet_node_numbering[1][0] = 1;
            facet_node_numbering[1][1] = 2;
            facet_node_numbering[2][0] = 2;
            facet_node_numbering[2][1] = 0;
        }
        else {
            throw_Exception(elements.name() + " is not \"Pentagon\", \"Quadrilateral\", or \"Triangle\"", Here());
        }

        std::vector<idx_t> facet_nodes(nb_nodes_in_facet);

        const idx_t e_start = range[t].start();
        const idx_t e_end   = range[t].end();

        for (idx_t e = e_start; e < e_end; ++e) {
            if (patch(e)) {
                continue;
            }

            for (idx_t f = 0; f < nb_facets_in_elem; ++f) {
                bool found_face = false;

                for (idx_t jnode = 0; jnode < nb_nodes_in_facet; ++jnode) {
                    facet_nodes[jnode] = elem_nodes(e, facet_node_numbering[f][jnode]);
                }

                int node = facet_nodes[0];
                for (idx_t face : node_to_facet[node]) {
                    idx_t nb_matched_nodes = 0;
                    if (nb_nodes_in_facet > 1)  // 2D or 3D
                    {
                        for (idx_t jnode = 0; jnode < nb_nodes_in_facet; ++jnode) {
                            idx_t other_node = facet_nodes[jnode];
                            for (idx_t other_face : node_to_facet[other_node]) {
                                if (other_face == face) {
                                    ++nb_matched_nodes;
                                    break;
                                }
                            }
                        }
                        if (nb_matched_nodes == nb_nodes_in_facet) {
                            connectivity_facet_to_elem[2 * face + 1] = e + elements.begin();
                            ++nb_inner_facets;
                            found_face = true;
                            break;
                        }
                    }
                }

                if (found_face == false) {
                    connectivity_facet_to_elem.emplace_back(elements.begin() + e);
                    // if 2nd element stays missing_value, it is a bdry face
                    connectivity_facet_to_elem.emplace_back(missing_value);
                    for (idx_t n = 0; n < nb_nodes_in_facet; ++n) {
                        node_to_facet[facet_nodes[n]].emplace_back(nb_facets);
                        facet_nodes_data.emplace_back(facet_nodes[n]);
                    }
                    ++nb_facets;
                }
            }
        }
    }
}

void accumulate_facets_ordered_by_halo(const mesh::HybridElements& cells, const mesh::Nodes& nodes,
                                       std::vector<idx_t>& facet_nodes_data,  // shape(nb_facets,nb_nodes_per_facet)
                                       std::vector<idx_t>& connectivity_facet_to_elem, idx_t& nb_facets,
                                       idx_t& nb_inner_facets, idx_t& missing_value, std::vector<idx_t>& halo_offsets) {
    ATLAS_TRACE();

    static int MAXHALO = 50;
    std::vector<std::vector<array::Range>> ranges(MAXHALO, std::vector<array::Range>(cells.nb_types()));

    int maxhalo{0};
    for (idx_t t = 0; t < cells.nb_types(); ++t) {
        const mesh::Elements& elements = cells.elements(t);
        auto elem_halo                 = elements.view<int, 1>(elements.halo());
        idx_t nb_elems                 = elements.size();

        int halo{0};
        int begin{0};
        int end{0};
        for (idx_t e = 0; e < nb_elems; ++e) {
            ATLAS_ASSERT(elem_halo(e) >= halo);
            if (elem_halo(e) > halo) {
                end             = e;
                ranges[halo][t] = array::Range{begin, end};
                begin           = end;
                ++halo;
            }
        }
        end             = nb_elems;
        ranges[halo][t] = array::Range{begin, end};
        maxhalo         = std::max(halo, maxhalo);
    }


    missing_value = -1;
    std::vector<std::vector<idx_t>> node_to_facet(nodes.size());
    for (auto& facets : node_to_facet) {
        facets.reserve(6);
    }
    nb_facets       = 0;
    nb_inner_facets = 0;


    halo_offsets = std::vector<idx_t>{0};
    for (int h = 0; h <= maxhalo; ++h) {
        accumulate_facets_in_range(ranges[h], cells, nodes, facet_nodes_data, connectivity_facet_to_elem, nb_facets,
                                   nb_inner_facets, missing_value, node_to_facet);
        halo_offsets.emplace_back(nb_facets);
    }
}


}  // namespace detail
}  // namespace mesh
}  // namespace atlas
