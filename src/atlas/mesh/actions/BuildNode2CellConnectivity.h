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

namespace atlas {
class Mesh;
}  // namespace atlas

namespace atlas {
namespace mesh {
namespace actions {

class BuildNode2CellConnectivity {
public:
    BuildNode2CellConnectivity(Mesh& mesh);
    void operator()();

private:
    Mesh& mesh_;
};

inline void build_node_to_cell_connectivity(Mesh& mesh) {
    BuildNode2CellConnectivity{mesh}();
}

}  // namespace actions
}  // namespace mesh
}  // namespace atlas

// ------------------------------------------------------------------

namespace atlas {
namespace mesh {
namespace detail {
class MeshImpl;
}
}  // namespace mesh
}  // namespace atlas

// C wrapper interfaces to C++ routines
extern "C" {
void atlas__build_node_to_cell_connectivity(atlas::mesh::detail::MeshImpl* mesh);
}

// ------------------------------------------------------------------
