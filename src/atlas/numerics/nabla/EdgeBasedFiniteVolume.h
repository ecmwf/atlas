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
    class EdgeBasedFiniteVolume;
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

  atlas::functionspace::EdgeBasedFiniteVolume const *fvm_;
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

class EdgeBasedFiniteVolume : public NodesFunctionSpace {
public:
  EdgeBasedFiniteVolume(Mesh &, const Halo & = Halo(1) );

  virtual std::string name() const { return "EdgeBasedFiniteVolume"; }

  const NodesFunctionSpace& nodes_fs() const { return *this; }
        NodesFunctionSpace& nodes_fs()       { return *this; }

private: // data

    atlas::FunctionSpace* edges_; // non-const because functionspace may modify mesh
};

}
}












#endif // atlas_numerics_nabla_EdgeBasedFiniteVolume_h
