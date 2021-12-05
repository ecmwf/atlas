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

#include <string>
#include <vector>

#include "atlas/library/config.h"

namespace atlas {
class Mesh;

namespace mesh {
namespace actions {

class BuildHalo {
public:
    BuildHalo(Mesh& mesh);
    void operator()(int nb_elems);

private:
    Mesh& mesh_;

public:
    std::vector<idx_t> periodic_points_local_index_;
    std::vector<std::vector<idx_t>> periodic_cells_local_index_;
};

/// @brief Enlarge each partition of the mesh with a halo of elements
/// @param [inout] mesh      The mesh to enlarge
/// @param [in]    nb_elems  Size of the halo
/// @author Willem Deconinck
/// @date June 2014
inline void build_halo(Mesh& mesh, int nb_elems) {
    BuildHalo f(mesh);
    f(nb_elems);
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" {
void atlas__build_halo(Mesh::Implementation* mesh, int nb_elems);
}
// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
