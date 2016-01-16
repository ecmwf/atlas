/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_AccumulateFaces_h
#define atlas_util_AccumulateFaces_h

#include "atlas/FunctionSpace.h"
#include "atlas/Parameters.h"

namespace atlas { class FunctionSpace; };
namespace atlas { namespace mesh { class HybridElements; } }
namespace atlas { namespace mesh { class Nodes; } }

namespace atlas {
namespace util {

struct Face
{
	ElementRef& operator[](size_t i) { return elems[i]; }
	bool is_bdry() const { return (elems[1].f < 0); }
	ElementRef elems[2];
};

void accumulate_faces(
		FunctionSpace& func_space,
		std::vector< std::vector<int> >& node_to_face,
		std::vector<int>& face_nodes_data,
		std::vector< Face >& connectivity_edge_to_elem,
		size_t& nb_faces,
		size_t& nb_inner_faces );

// currently only supports 2D meshes. Little work needed for 3D.
void accumulate_facets(
    const mesh::HybridElements &cells,
    const mesh::Nodes &nodes,
    std::vector< idx_t > &facet_nodes_data, // shape(nb_facets,nb_nodes_per_facet)
    std::vector< idx_t > &connectivity_facet_to_elem,
    size_t &nb_facets,
    size_t &nb_inner_facets );

} // namespace util
} // namespace atlas

#endif
