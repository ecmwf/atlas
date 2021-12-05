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

#include <vector>

#include "atlas/library/config.h"
#include "atlas/numerics/Nabla.h"

namespace atlas {
namespace numerics {
namespace fvm {
class Method;
}
}  // namespace numerics
}  // namespace atlas

namespace atlas {
class Field;
}

namespace atlas {
namespace numerics {
namespace fvm {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Nabla : public atlas::numerics::NablaImpl {
public:
    Nabla(const atlas::numerics::Method&, const eckit::Parametrisation&);
    virtual ~Nabla() override;

    virtual void gradient(const Field& scalar, Field& grad) const override;
    virtual void divergence(const Field& vector, Field& div) const override;
    virtual void curl(const Field& vector, Field& curl) const override;
    virtual void laplacian(const Field& scalar, Field& laplacian) const override;

    virtual const FunctionSpace& functionspace() const override;

private:
    void setup();

    void gradient_of_scalar(const Field& scalar, Field& grad) const;
    void gradient_of_vector(const Field& vector, Field& grad) const;

private:
    fvm::Method const* fvm_;
    std::vector<idx_t> pole_edges_;
    int metric_approach_{0};
};
#endif
// ------------------------------------------------------------------

}  // namespace fvm
}  // namespace numerics
}  // namespace atlas
