/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <vector>

#include "atlas/library/config.h"

namespace atlas {
namespace mesh {
class HybridElements;
}
}  // namespace atlas
namespace atlas {
namespace mesh {
class Nodes;
}
}  // namespace atlas

namespace atlas {
namespace mesh {
namespace detail {

// currently only supports 2D meshes. Little work needed for 3D.
void accumulate_facets(const mesh::HybridElements& cells, const mesh::Nodes& nodes,
                       std::vector<idx_t>& facet_nodes_data,  // shape(nb_facets,nb_nodes_per_facet)
                       std::vector<idx_t>& connectivity_facet_to_elem, idx_t& nb_facets, idx_t& nb_inner_facets,
                       idx_t& missing_value);

// currently only supports 2D meshes. Little work needed for 3D.
void accumulate_facets_ordered_by_halo(const mesh::HybridElements& cells, const mesh::Nodes& nodes,
                                       std::vector<idx_t>& facet_nodes_data,  // shape(nb_facets,nb_nodes_per_facet)
                                       std::vector<idx_t>& connectivity_facet_to_elem, idx_t& nb_facets,
                                       idx_t& nb_inner_facets, idx_t& missing_value, std::vector<idx_t>& halo_offsets);

}  // namespace detail
}  // namespace mesh
}  // namespace atlas
