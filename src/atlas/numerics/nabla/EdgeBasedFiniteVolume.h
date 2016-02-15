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

#include <vector>
#include "atlas/numerics/Nabla.h"
#include "eckit/memory/SharedPtr.h"

namespace atlas {  namespace functionspace {  class EdgeBasedFiniteVolume; } }
namespace atlas {  class Field; }

namespace atlas {
namespace numerics {
namespace nabla {

class EdgeBasedFiniteVolume : public Nabla {

public:
  EdgeBasedFiniteVolume(const FunctionSpace &, const eckit::Parametrisation &);
  virtual ~EdgeBasedFiniteVolume();

  void gradient(const Field &scalar, Field &grad) const;
  void divergence(const Field &vector, Field &div) const;
  void curl(const Field &vector, Field &curl) const;
  void laplacian(const Field &scalar, Field &laplacian) const;

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

#endif // atlas_numerics_nabla_EdgeBasedFiniteVolume_h
