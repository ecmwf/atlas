/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_numerics_nabla_EdgeBasedFiniteVolume_h
#define atlas_numerics_nabla_EdgeBasedFiniteVolume_h

#include "atlas/numerics/Nabla.h"
#include "eckit/memory/SharedPtr.h"

namespace atlas {
  namespace functionspace {
    class NodesFunctionSpace;
    class EdgeBasedFiniteVolumeFunctionSpace;
  }
  class Field;
}

namespace atlas {
namespace numerics {
namespace nabla {

class EdgeBasedFiniteVolume : public Nabla {

public:
  EdgeBasedFiniteVolume(const next::FunctionSpace &, const eckit::Parametrisation &);
  virtual ~EdgeBasedFiniteVolume();

  void gradient(const Field &field, Field &grad);

private:
  void setup();

private:

  atlas::functionspace::EdgeBasedFiniteVolumeFunctionSpace const *fvm_;
  std::vector<size_t> pole_edges_;
};

// ------------------------------------------------------------------

} // namespace nabla
} // namespace numerics
} // namespace atlas





#include "atlas/FunctionSpace.h"
#include "atlas/functionspace/NodesFunctionSpace.h"


namespace atlas {
class Mesh;
namespace functionspace {
class EdgeBasedFiniteVolumeFunctionSpace : public next::FunctionSpace {
  friend class atlas::numerics::nabla::EdgeBasedFiniteVolume;
public:
  EdgeBasedFiniteVolumeFunctionSpace(Mesh &, const Halo & = Halo(1) );

private: // data

    Mesh& mesh_; // non-const because functionspace may modify mesh
    eckit::SharedPtr<atlas::functionspace::NodesColumnFunctionSpace> nodes_; // non-const because functionspace may modify mesh
    atlas::FunctionSpace* edges_; // non-const because functionspace may modify mesh
    Halo halo_;
    size_t nb_edges_;
    size_t nb_edges_global_;
    std::vector<size_t> nb_edges_global_foreach_rank_;
};
}
}












#endif // atlas_numerics_nabla_EdgeBasedFiniteVolume_h
