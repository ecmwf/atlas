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
#include "atlas/functionspace/Nodes.h"
#include "atlas/functionspace/Edges.h"
#include "atlas/numerics/Method.h"

namespace eckit { class Parametrisation; }
namespace atlas { class Mesh; }
namespace atlas { namespace mesh { class HybridElements; } }
namespace atlas { namespace mesh { class Nodes; } }

namespace atlas {
namespace numerics {
namespace fvm {

class Method : public numerics::Method {

public:

  Method(Mesh &, const eckit::Parametrisation &);
  Method(Mesh &, const mesh::Halo &);
  Method(Mesh &);

  virtual std::string name() const { return "fvm"; }

  const atlas::Mesh& mesh() const { return mesh_; }
        atlas::Mesh& mesh()       { return mesh_; }

  const mesh::Nodes& nodes() const { return nodes_; }
        mesh::Nodes& nodes()       { return nodes_; }

  const mesh::HybridElements& edges() const { return edges_; }
        mesh::HybridElements& edges()       { return edges_; }

  const functionspace::Nodes& nodes_fs() const { return *nodes_fs_; }
        functionspace::Nodes& nodes_fs()       { return *nodes_fs_; }

  const functionspace::Edges& edges_fs() const { return *edges_fs_; }
        functionspace::Edges& edges_fs()       { return *edges_fs_; }

  const double& radius() const { return radius_; }

private:

  void setup();

private: // data

    mesh::Halo             halo_;
    Mesh                   &mesh_; // non-const because functionspace may modify mesh
    mesh::Nodes            &nodes_;
    mesh::HybridElements   &edges_;
    eckit::SharedPtr<functionspace::Nodes>    nodes_fs_;
    eckit::SharedPtr<functionspace::Edges>    edges_fs_;

    double radius_;


};


// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
#define Parametrisation eckit::Parametrisation
#define functionspace_Nodes functionspace::Nodes
extern "C"
{
  Method* atlas__numerics__fvm__Method__new (Mesh* mesh, const Parametrisation* params);
  functionspace_Nodes* atlas__numerics__fvm__Method__nodes_fs (Method* This);
}
#undef Parametrisation
#undef functionspace_Nodes

} // namespace fvm
} // namespace numerics
} // namespace atlas


#endif // atlas_numerics_fvm_Method_h
