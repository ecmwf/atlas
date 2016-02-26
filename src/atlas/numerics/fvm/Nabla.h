/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_numerics_fvm_Nabla_h
#define atlas_numerics_fvm_Nabla_h

#include <vector>
#include "eckit/memory/SharedPtr.h"
#include "atlas/numerics/Nabla.h"

namespace atlas {
namespace numerics {
namespace fvm  {
    class Method;
} } }

namespace atlas {namespace field { class Field; } }

namespace atlas {
namespace numerics {
namespace fvm {

class Nabla : public atlas::numerics::Nabla {

public:
  Nabla(const atlas::numerics::Method &, const eckit::Parametrisation &);
  virtual ~Nabla();

  void gradient(const field::Field &scalar, field::Field &grad) const;
  void divergence(const field::Field &vector, field::Field &div) const;
  void curl(const field::Field &vector, field::Field &curl) const;
  void laplacian(const field::Field &scalar, field::Field &laplacian) const;

private:
  void setup();

private:

  fvm::Method const *fvm_;
  std::vector<size_t> pole_edges_;
};

// ------------------------------------------------------------------

} // namespace fvm
} // namespace numerics
} // namespace atlas

#endif // atlas_numerics_fvm_Nabla_h
