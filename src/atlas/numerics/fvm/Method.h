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

#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/numerics/Method.h"

namespace eckit {
class Parametrisation;
}
namespace atlas {
class Mesh;
}
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
namespace numerics {
namespace fvm {

class Method : public numerics::Method {
public:
    Method(Mesh&, const eckit::Configuration&);
    Method(Mesh&, const mesh::Halo&);
    Method(Mesh&);

    virtual std::string& name() const override {
        static std::string _name{"fvm"};
        return _name;
    }

    const atlas::Mesh& mesh() const { return mesh_; }
    atlas::Mesh& mesh() { return mesh_; }

    size_t levels() const { return levels_; }

    const functionspace::NodeColumns& node_columns() const { return node_columns_; }
    functionspace::NodeColumns& node_columns() { return node_columns_; }

    const functionspace::EdgeColumns& edge_columns() const { return edge_columns_; }
    functionspace::EdgeColumns& edge_columns() { return edge_columns_; }

    const double& radius() const { return radius_; }

private:
    void setup();

private:         // data
    Mesh mesh_;  // non-const because functionspace may modify mesh
    size_t levels_;
    mesh::Halo halo_;
    mesh::Nodes& nodes_;
    mesh::HybridElements& edges_;
    functionspace::NodeColumns node_columns_;
    functionspace::EdgeColumns edge_columns_;

    double radius_;
};

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines

extern "C" {
Method* atlas__numerics__fvm__Method__new(Mesh::Implementation* mesh, const eckit::Configuration* params);
const functionspace::detail::NodeColumns* atlas__numerics__fvm__Method__functionspace_nodes(Method* This);
const functionspace::detail::EdgeColumns* atlas__numerics__fvm__Method__functionspace_edges(Method* This);
}

}  // namespace fvm
}  // namespace numerics
}  // namespace atlas
