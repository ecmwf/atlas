/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_numerics_fvm_Method_h
#define atlas_numerics_fvm_Method_h

#include <string>
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/EdgeColumns.h"
#include "atlas/numerics/Method.h"

namespace eckit { class Parametrisation; }
namespace atlas { namespace mesh { class Mesh; } }
namespace atlas { namespace mesh { class HybridElements; } }
namespace atlas { namespace mesh { class Nodes; } }

namespace atlas {
namespace numerics {
namespace fvm {

class Method : public numerics::Method {

public:

  Method(mesh::Mesh &, const eckit::Parametrisation &);
  Method(mesh::Mesh &, const mesh::Halo &);
  Method(mesh::Mesh &);

  virtual std::string name() const { return "fvm"; }

  const atlas::mesh::Mesh& mesh() const { return mesh_; }
        atlas::mesh::Mesh& mesh()       { return mesh_; }

  const functionspace::NodeColumns& node_columns() const { return *node_columns_; }
        functionspace::NodeColumns& node_columns()       { return *node_columns_; }

  const functionspace::EdgeColumns& edge_columns() const { return *edge_columns_; }
        functionspace::EdgeColumns& edge_columns()       { return *edge_columns_; }

  const double& radius() const { return radius_; }

private:

  void setup();

private: // data

    mesh::Mesh             &mesh_; // non-const because functionspace may modify mesh
    mesh::Halo             halo_;
    mesh::Nodes            &nodes_;
    mesh::HybridElements   &edges_;
    eckit::SharedPtr<functionspace::NodeColumns>  node_columns_;
    eckit::SharedPtr<functionspace::EdgeColumns>  edge_columns_;

    double radius_;


};


// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
#define Parametrisation eckit::Parametrisation
#define mesh_Mesh mesh::Mesh
#define functionspace_NodeColumns functionspace::NodeColumns
#define functionspace_EdgeColumns functionspace::EdgeColumns

extern "C"
{
  Method* atlas__numerics__fvm__Method__new (mesh_Mesh* mesh, const Parametrisation* params);
  functionspace_NodeColumns* atlas__numerics__fvm__Method__functionspace_nodes (Method* This);
  functionspace_EdgeColumns* atlas__numerics__fvm__Method__functionspace_edges (Method* This);
}
#undef Parametrisation
#undef mesh_Mesh
#undef functionspace_NodeColumns
#undef functionspace_EdgeColumns

} // namespace fvm
} // namespace numerics
} // namespace atlas


#endif // atlas_numerics_fvm_Method_h
